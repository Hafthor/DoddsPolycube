namespace DoddsPolycube;

using Num = ulong;

// SIMD Vectorization Analysis for CountExtensionsSubsetUnsafe2
// ==============================================================
//
// The hot loop accesses 6 neighbors at byte offsets ±X(65), ±Y(66), ±Z(80) — scattered,
// non-contiguous locations in a 1440-byte board that fits entirely in L1 cache.
//
// Why classical SIMD doesn't help:
// - ARM NEON (M2 Max): No hardware byte gather/scatter. Loading 6 scattered bytes into a
//   vector requires 6 individual loads + inserts, which is slower than just using them scalar.
// - AVX2 (x86): Has VPGATHERDD but only for 32/64-bit elements, not bytes.
// - The 6 neighbor values are at DIFFERENT non-uniform offsets, so no contiguous vector load.
// - The board is only 1440 bytes → already in L1 cache → no prefetch benefit.
// - Branches in the >127 checks are well-predicted → branchless (mask*value) is slower.
// - The M2's out-of-order engine already extracts maximum ILP from the scalar code.
//
// Approaches tested and their results (filter 0, N=16, M2 Max):
// - Baseline (CountExtensionsSubsetUnsafe2):           ~6.3s
// - Branchless with mask multiply:                     ~7.2s (worse - extra multiply overhead)
// - Speculative load all 30 sub-neighbors + masks:     ~7.2s (worse - wasted loads)
// - Template-specialized recursion per depth:          ~7.7s (worse - closure overhead)
// - 2x unrolled Work4 loop:                            ~6.3s (same - OoO already extracts ILP)
// - Software prefetch in Work4:                        ~6.5s (same - data already in L1)
// - Pre-loaded 6 neighbor values before branching:     ~6.3s (same - compiler already does this)
// - Cached pointer addresses for Wind/Unwind:          ~6.4s (same - CSE already handles this)
//
// The only path to significant speedup is ALGORITHMIC:
// - Extend the depth==4 analytical base case to depth==5 (eliminate one recursion level)
//   This would reduce total Work4 invocations by the branching factor (~5-8x)
//   But requires deriving a 4th-order combinatorial correction formula.
// - NativeAOT compilation may yield 5-15% from better code generation.
//
// The two methods below are kept for reference as the best micro-optimization attempts.

public partial class Program {
    /// <summary>
    /// Pre-loads all 6 neighbor values before any branching in Work4.
    /// Theory: issuing all 6 byte loads before the first branch lets the CPU's load unit
    /// dispatch them in parallel. In practice, the JIT and OoO engine already achieve this.
    /// Performance: within noise of baseline (~6.3s).
    /// </summary>
    private static unsafe Num CountExtensionsSubsetBranchPreload(int filter) {
        if (quiting) return 0;
        byte* byteBoard = stackalloc byte[(N + 2) * Z];
        byte** refStack = stackalloc byte*[(N - 2) * 4];
        *refStack = byteBoard += Z;

        for (byte* i = byteBoard + (N + 1) * Z; --i != byteBoard;)
            *i = 255;

        return CountExtensions(N, refStack + 1, refStack + (N - 2) * 4);

        Num CountExtensions(int depth, byte** stackTop1, byte** stackTop2) {
            Num count = 0;
            byte** stackTopOriginal = stackTop1;
            while (stackTop1 != refStack) {
                byte* index = *--stackTop1;
                if (depth != FilterDepth || stackTop1 - refStack == filter) {
                    byte** stackTopInner = stackTop1;

                    if (++*(index - Z) == 0) *stackTopInner++ = index - Z;
                    if (++*(index - Y) == 0) *stackTopInner++ = index - Y;
                    if (++*(index - X) == 0) *stackTopInner++ = index - X;
                    if (++*(index + X) == 0) *stackTopInner++ = index + X;
                    if (++*(index + Y) == 0) *stackTopInner++ = index + Y;
                    if (++*(index + Z) == 0) *stackTopInner++ = index + Z;

                    if (depth == 4) {
                        byte** stackTop = stackTopInner, stackTopTemp = stackTopInner;
                        int length = (int)(stackTop - refStack), lengthPlus = (length << 1) - 511;
                        count += (Num)(length * (length - 1) * (length - 2) / 6);
                        for (; stackTopTemp != refStack;) {
                            byte* i = *--stackTopTemp;
                            // Pre-load all 6 neighbor values to issue loads early
                            byte vZ0 = *(i - Z), vY0 = *(i - Y), vX0 = *(i - X);
                            byte vX1 = *(i + X), vY1 = *(i + Y), vZ1 = *(i + Z);
                            int neighbours = 0, subCount = 128, localCount = 0;
                            if (vZ0 > 127) {
                                *(i - Z) = --vZ0;
                                localCount += vZ0;
                                subCount += *(i - Z - Z) + *(i - Z - X) + *(i - Z - Y) + *(i - Z + X) + *(i - Z + Y);
                                neighbours++;
                            }
                            if (vY0 > 127) {
                                *(i - Y) = --vY0;
                                localCount += vY0;
                                subCount += *(i - Y - Y) + *(i - Y - X) + *(i - Y - Z) + *(i - Y + X) + *(i - Y + Z);
                                neighbours++;
                            }
                            if (vX0 > 127) {
                                *(i - X) = --vX0;
                                localCount += vX0;
                                subCount += *(i - X - X) + *(i - X - Y) + *(i - X - Z) + *(i - X + Y) + *(i - X + Z);
                                neighbours++;
                            }
                            if (vX1 > 127) {
                                *(i + X) = --vX1;
                                localCount += vX1;
                                subCount += *(i + X + X) + *(i + X + Y) + *(i + X + Z) + *(i + X - Y) + *(i + X - Z);
                                neighbours++;
                            }
                            if (vY1 > 127) {
                                *(i + Y) = --vY1;
                                localCount += vY1;
                                subCount += *(i + Y + Y) + *(i + Y + X) + *(i + Y + Z) + *(i + Y - X) + *(i + Y - Z);
                                neighbours++;
                            }
                            if (vZ1 > 127) {
                                *(i + Z) = --vZ1;
                                localCount += vZ1;
                                subCount += *(i + Z + Z) + *(i + Z + X) + *(i + Z + Y) + *(i + Z - X) + *(i + Z - Y);
                                neighbours++;
                            }
                            count += (Num)localCount +
                                     (Num)((subCount >> 8) + (neighbours * (neighbours + lengthPlus) >> 1));
                        }
                        while (stackTop != refStack) {
                            byte* i = *--stackTop;
                            *(i - Z) |= (byte)(*(i - Z) >> 4);
                            *(i - Y) |= (byte)(*(i - Y) >> 4);
                            *(i - X) |= (byte)(*(i - X) >> 4);
                            *(i + X) |= (byte)(*(i + X) >> 4);
                            *(i + Y) |= (byte)(*(i + Y) >> 4);
                            *(i + Z) |= (byte)(*(i + Z) >> 4);
                        }
                    } else {
                        count += CountExtensions(depth - 1, stackTopInner, stackTop2);
                    }

                    --*(index - Z);
                    --*(index - Y);
                    --*(index - X);
                    --*(index + X);
                    --*(index + Y);
                    --*(index + Z);
                }

                *--stackTop2 = index;
            }
            while (stackTop1 != stackTopOriginal)
                *stackTop1++ = *stackTop2++;
            return count;
        }
    }

    /// <summary>
    /// Caches the 6 neighbor pointer addresses computed during Wind for reuse in Unwind.
    /// Theory: avoids recomputing index±X/Y/Z in the Unwind phase. In practice, the JIT
    /// already performs common subexpression elimination (CSE) for these.
    /// Performance: within noise of baseline (~6.4s).
    /// </summary>
    private static unsafe Num CountExtensionsSubsetSimd2(int filter) {
        if (quiting) return 0;
        byte* byteBoard = stackalloc byte[(N + 2) * Z];
        byte** refStack = stackalloc byte*[(N - 2) * 4];
        *refStack = byteBoard += Z;

        for (byte* i = byteBoard + (N + 1) * Z; --i != byteBoard;)
            *i = 255;

        return CountExtensions(N, refStack + 1, refStack + (N - 2) * 4);

        Num CountExtensions(int depth, byte** stackTop1, byte** stackTop2) {
            Num count = 0;
            byte** stackTopOriginal = stackTop1;
            while (stackTop1 != refStack) {
                byte* index = *--stackTop1;
                if (depth != FilterDepth || stackTop1 - refStack == filter) {
                    byte** stackTopInner = stackTop1;

                    // Cache neighbor addresses for Wind + Unwind
                    byte* nZ0 = index - Z, nY0 = index - Y, nX0 = index - X;
                    byte* nX1 = index + X, nY1 = index + Y, nZ1 = index + Z;
                    if (++*nZ0 == 0) *stackTopInner++ = nZ0;
                    if (++*nY0 == 0) *stackTopInner++ = nY0;
                    if (++*nX0 == 0) *stackTopInner++ = nX0;
                    if (++*nX1 == 0) *stackTopInner++ = nX1;
                    if (++*nY1 == 0) *stackTopInner++ = nY1;
                    if (++*nZ1 == 0) *stackTopInner++ = nZ1;

                    if (depth == 4) {
                        byte** stackTop = stackTopInner, stackTopTemp = stackTopInner;
                        int length = (int)(stackTop - refStack), lengthPlus = (length << 1) - 511;
                        count += (Num)(length * (length - 1) * (length - 2) / 6);
                        for (; stackTopTemp != refStack;) {
                            byte* i = *--stackTopTemp;
                            int neighbours = 0, subCount = 128, localCount = 0;
                            if (*(i - Z) > 127) {
                                localCount += --*(i - Z);
                                subCount += *(i - Z - Z) + *(i - Z - X) + *(i - Z - Y) + *(i - Z + X) + *(i - Z + Y);
                                neighbours++;
                            }
                            if (*(i - Y) > 127) {
                                localCount += --*(i - Y);
                                subCount += *(i - Y - Y) + *(i - Y - X) + *(i - Y - Z) + *(i - Y + X) + *(i - Y + Z);
                                neighbours++;
                            }
                            if (*(i - X) > 127) {
                                localCount += --*(i - X);
                                subCount += *(i - X - X) + *(i - X - Y) + *(i - X - Z) + *(i - X + Y) + *(i - X + Z);
                                neighbours++;
                            }
                            if (*(i + X) > 127) {
                                localCount += --*(i + X);
                                subCount += *(i + X + X) + *(i + X + Y) + *(i + X + Z) + *(i + X - Y) + *(i + X - Z);
                                neighbours++;
                            }
                            if (*(i + Y) > 127) {
                                localCount += --*(i + Y);
                                subCount += *(i + Y + Y) + *(i + Y + X) + *(i + Y + Z) + *(i + Y - X) + *(i + Y - Z);
                                neighbours++;
                            }
                            if (*(i + Z) > 127) {
                                localCount += --*(i + Z);
                                subCount += *(i + Z + Z) + *(i + Z + X) + *(i + Z + Y) + *(i + Z - X) + *(i + Z - Y);
                                neighbours++;
                            }
                            count += (Num)localCount +
                                     (Num)((subCount >> 8) + (neighbours * (neighbours + lengthPlus) >> 1));
                        }
                        while (stackTop != refStack) {
                            byte* i = *--stackTop;
                            *(i - Z) |= (byte)(*(i - Z) >> 4);
                            *(i - Y) |= (byte)(*(i - Y) >> 4);
                            *(i - X) |= (byte)(*(i - X) >> 4);
                            *(i + X) |= (byte)(*(i + X) >> 4);
                            *(i + Y) |= (byte)(*(i + Y) >> 4);
                            *(i + Z) |= (byte)(*(i + Z) >> 4);
                        }
                    } else {
                        count += CountExtensions(depth - 1, stackTopInner, stackTop2);
                    }

                    // Unwind using cached pointers
                    --*nZ0;
                    --*nY0;
                    --*nX0;
                    --*nX1;
                    --*nY1;
                    --*nZ1;
                }

                *--stackTop2 = index;
            }
            while (stackTop1 != stackTopOriginal)
                *stackTop1++ = *stackTop2++;
            return count;
        }
    }
}





