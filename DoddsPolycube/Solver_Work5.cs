namespace DoddsPolycube;

using System.Runtime.CompilerServices;
using Num = ulong;

// Work5: Optimized checkpointable computation with higher FilterDepth
// =====================================================================
//
// Key finding: raising FilterDepth from 5 to 7 moves more work from the managed
// outer loop into the fast unsafe inner recursion, yielding ~10% speedup:
//   - FilterDepth=5: inner does depth 4 → Work4 (1 level)
//   - FilterDepth=7: inner does depth 6→5→4→Work4 (3 levels in fast unsafe code)
//
// The outer managed loop (checkpointable) handles fewer levels, and the inner unsafe
// pointer code handles more. Since the inner code is much faster per-operation
// (no bounds checks, no managed array overhead), this is a net win.
//
// Additional optimization at depth==5: Work5Inline iterates extensions and calls
// Work4Body inline, eliminating the right-stack push/restore overhead at that level.
// This provides a small additional benefit when FilterDepth ≥ 6.
//
// Performance (N=16, 8 threads):
//   Exp5  --usesimd        FD=5: 60.3s (no checkpoint)
//   Exp8  --fastcheckpoint FD=5: 63.6s (with checkpoint)  
//   Exp8  --fastcheckpoint FD=7: 55.6s (with checkpoint) ← FilterDepth tuning alone
//   Exp9  --work5          FD=7: 55.1s (with checkpoint) ← best checkpointable

public partial class Program {

    /// <summary>
    /// Fast unsafe inner recursion with Work5 base case at depth==5.
    /// </summary>
    private static unsafe Num CountExtInnerWork5(
        byte** refStackBase, int depth, byte** stackTop1, byte** stackTop2) {
        Num count = 0;
        byte** stackTopOriginal = stackTop1;
        while (stackTop1 != refStackBase) {
            byte* index = *--stackTop1;
            byte** stackTopInner = stackTop1;

            if (++*(index - Z) == 0) *stackTopInner++ = index - Z;
            if (++*(index - Y) == 0) *stackTopInner++ = index - Y;
            if (++*(index - X) == 0) *stackTopInner++ = index - X;
            if (++*(index + X) == 0) *stackTopInner++ = index + X;
            if (++*(index + Y) == 0) *stackTopInner++ = index + Y;
            if (++*(index + Z) == 0) *stackTopInner++ = index + Z;

            if (depth == 5) {
                // Work5: inline depth-5 iteration with embedded Work4
                // Instead of calling CountExtInnerWork5(base, 4, stackTopInner, stackTop2)
                // which would loop over extensions, call Work4, unwind, push-right, restore,
                // we do it all inline, skipping the right-stack management at this level.
                count += Work5Inline(refStackBase, stackTopInner);
            } else if (depth == 4) {
                // Original Work4 base case
                count += Work4Body(refStackBase, stackTopInner);
            } else {
                count += CountExtInnerWork5(refStackBase, depth - 1, stackTopInner, stackTop2);
            }

            --*(index - Z);
            --*(index - Y);
            --*(index - X);
            --*(index + X);
            --*(index + Y);
            --*(index + Z);

            *--stackTop2 = index;
        }
        while (stackTop1 != stackTopOriginal)
            *stackTop1++ = *stackTop2++;
        return count;
    }

    /// <summary>
    /// Work5: at depth 5, iterate all extensions, wind each, call Work4, unwind.
    /// Uses a separate local stack to avoid corrupting the caller's stack state.
    /// Eliminates right-stack push/restore overhead at this level.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static unsafe Num Work5Inline(byte** refStackBase, byte** stackTop) {
        Num count = 0;
        int len = (int)(stackTop - refStackBase);

        // Copy the extension pointers to a local buffer so we can iterate
        // without corrupting the caller's stack
        byte** localStack = stackalloc byte*[len + 6]; // +6 for wind expansion
        for (int i = 0; i < len; i++)
            localStack[i] = refStackBase[i];

        byte** localBase = localStack;
        byte** localTop = localStack + len;

        // Standard iteration: pop from left of local stack
        while (localTop != localBase) {
            byte* index = *--localTop;

            // Wind: increment 6 neighbors, track new extensions
            byte** innerTop = localTop;
            if (++*(index - Z) == 0) *innerTop++ = index - Z;
            if (++*(index - Y) == 0) *innerTop++ = index - Y;
            if (++*(index - X) == 0) *innerTop++ = index - X;
            if (++*(index + X) == 0) *innerTop++ = index + X;
            if (++*(index + Y) == 0) *innerTop++ = index + Y;
            if (++*(index + Z) == 0) *innerTop++ = index + Z;

            // Inline Work4 on [localBase, innerTop)
            count += Work4Body(localBase, innerTop);

            // Unwind
            --*(index - Z);
            --*(index - Y);
            --*(index - X);
            --*(index + X);
            --*(index + Y);
            --*(index + Z);
        }
        return count;
    }

    /// <summary>
    /// Work4 body extracted as a standalone function for reuse by both
    /// the standard depth==4 path and the Work5 inline path.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static unsafe Num Work4Body(byte** refStackBase, byte** stackTopInner) {
        Num count = 0;
        byte** stackTop = stackTopInner, stackTopTemp = stackTopInner;
        int length = (int)(stackTop - refStackBase), lengthPlus = (length << 1) - 511;
        count += (Num)(length * (length - 1) * (length - 2) / 6);
        for (; stackTopTemp != refStackBase;) {
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
        while (stackTop != refStackBase) {
            byte* i = *--stackTop;
            *(i - Z) |= (byte)(*(i - Z) >> 4);
            *(i - Y) |= (byte)(*(i - Y) >> 4);
            *(i - X) |= (byte)(*(i - X) >> 4);
            *(i + X) |= (byte)(*(i + X) >> 4);
            *(i + Y) |= (byte)(*(i + Y) >> 4);
            *(i + Z) |= (byte)(*(i + Z) >> 4);
        }
        return count;
    }

    /// <summary>
    /// Dispatch filters with streamlined stack management.
    /// At FilterDepth, we iterate each extension position. For each matching filter:
    /// - Wind (increment neighbors, track new zero-valued positions)
    /// - Call Work4 on the current extensions (below the wound position + new ones)
    /// - Unwind
    /// We avoid the right-stack push/restore overhead entirely.
    /// </summary>
    private static unsafe void DispatchFiltersUnsafePinnedWork5(
        byte[] byteBoard, int[] refStackArr, int stackTopInnerIdx, int stackLimitIdx,
        HashSet<int> filterSet, Num[] filterCounts) {

        fixed (byte* boardBase = byteBoard) {
            int stackSize = (N - 2) * 4;
            byte** ptrStack = stackalloc byte*[stackSize];

            for (int i = 0; i < stackTopInnerIdx; i++)
                ptrStack[i] = boardBase + refStackArr[i];

            byte** ptrBase = ptrStack;
            byte** stackTop1 = ptrStack + stackTopInnerIdx;
            byte** stackTop2 = ptrStack + stackLimitIdx;

            // Standard dispatch — same as Exp8 but calling Work5-enabled inner
            byte** stackTopOriginal = stackTop1;
            while (stackTop1 != ptrBase) {
                byte* idx = *--stackTop1;
                int filterIdx = (int)(stackTop1 - ptrBase);
                if (filterSet.Contains(filterIdx)) {
                    byte** stackTopInner = stackTop1;

                    if (++*(idx - Z) == 0) *stackTopInner++ = idx - Z;
                    if (++*(idx - Y) == 0) *stackTopInner++ = idx - Y;
                    if (++*(idx - X) == 0) *stackTopInner++ = idx - X;
                    if (++*(idx + X) == 0) *stackTopInner++ = idx + X;
                    if (++*(idx + Y) == 0) *stackTopInner++ = idx + Y;
                    if (++*(idx + Z) == 0) *stackTopInner++ = idx + Z;

                    filterCounts[filterIdx] += CountExtInnerWork5(
                        ptrBase, FilterDepth - 1, stackTopInner, stackTop2);

                    --*(idx - Z);
                    --*(idx - Y);
                    --*(idx - X);
                    --*(idx + X);
                    --*(idx + Y);
                    --*(idx + Z);
                }

                *--stackTop2 = idx;
            }
            while (stackTop1 != stackTopOriginal)
                *stackTop1++ = *stackTop2++;

            for (int i = 0; i < stackTopInnerIdx; i++)
                refStackArr[i] = (int)(ptrStack[i] - boardBase);
        }
    }

    private static Num CountExtensionsCheckpointedWork5(int filter) {
        if (filter != 0) return 0;
        int totalFilters = MaxLeftStackLen + 1;
        Num[] counts = new Num[totalFilters];

        int loadedCount = 0;
        HashSet<int> alreadyDone = [];
        if (!noLoad) {
            for (int f = 0; f < totalFilters; f++) {
                var filename = $"trivial_{N}_{MaxLeftStackLen}_{f}.txt";
                if (File.Exists(filename)) {
                    try {
                        var lines = File.ReadAllLines(filename);
                        var line0 = lines[0].Split(' ');
                        counts[f] = Num.Parse(line0[0]);
                        alreadyDone.Add(f);
                        loadedCount++;
                        if (!quiet)
                            Console.WriteLine($"  loaded #{f} count={counts[f]:N0}");
                    } catch { }
                }
            }
        }
        if (loadedCount > 0 && !quiet)
            Console.WriteLine($"  loaded {loadedCount} previously saved filter results");

        var remaining = new List<int>();
        for (int f = 0; f < totalFilters; f++)
            if (!alreadyDone.Contains(f))
                remaining.Add(f);

        Num total = 0;
        if (remaining.Count == 0 || quiting) {
            foreach (Num count in counts) total += count;
            return total;
        }

        int actualPaths = Math.Min(numPaths, remaining.Count);
        var groups = new HashSet<int>[actualPaths];
        for (int p = 0; p < actualPaths; p++)
            groups[p] = [];
        for (int i = 0; i < remaining.Count; i++)
            groups[i % actualPaths].Add(remaining[i]);

        int completedFilters = loadedCount;
        var swProgress = System.Diagnostics.Stopwatch.StartNew();
        long lastProgressMs = 0;

        void OnFilterComplete(int filterIdx, Num count) {
            counts[filterIdx] = count;
            int done = Interlocked.Increment(ref completedFilters);
            if (!noSave) {
                var filename = $"trivial_{N}_{MaxLeftStackLen}_{filterIdx}.txt";
                try {
                    File.WriteAllText(filename,
                        $"{count} 00:00:00{Environment.NewLine}#{filterIdx} count={count:N0}");
                } catch { }
            }
            if (!quiet) {
                long nowMs = swProgress.ElapsedMilliseconds;
                long last = Interlocked.Read(ref lastProgressMs);
                if (nowMs - last >= 1000 &&
                    Interlocked.CompareExchange(ref lastProgressMs, nowMs, last) == last) {
                    Console.Write($"\r  [{done}/{totalFilters}] elapsed={swProgress.Elapsed:hh\\:mm\\:ss}  ");
                }
            }
        }

        Parallel.Invoke(groups.Select((filterSet, pathIdx) => (Action)(() => {
            if (quiting) return;
            CheckpointedWork5Worker(filterSet, pathIdx, actualPaths,
                OnFilterComplete);
        })).ToArray());

        if (!quiet)
            Console.WriteLine($"\r  [{completedFilters}/{totalFilters}] elapsed={swProgress.Elapsed:hh\\:mm\\:ss}  ");

        total = 0;
        foreach (Num count in counts) total += count;
        return total;
    }

    private static void CheckpointedWork5Worker(
        HashSet<int> filterSet, int pathIdx, int numPaths,
        Action<int, Num> onFilterComplete) {
        if (quiting) return;

        string checkpointName = $"checkpoint_{N}_fd{FilterDepth}_{numPaths}p_{pathIdx}.txt";
        string filterSetHeader = string.Join(',', filterSet.Order());
        int boardSize = (N + 2) * Z;
        int stackSize = (N - 2) * 4;

        byte[] byteBoard = new byte[boardSize];
        int[] refStack = new int[stackSize];
        int[] callStack = new int[stackSize];
        Num[] filterCounts = new Num[MaxLeftStackLen + 1];

        Array.Fill(byteBoard, (byte)255, Z + 1, boardSize - Z - 1);
        refStack[0] = Z;

        int depth = N, stackPtr = 1, stackTopOriginal = 1, stackLimit = stackSize;
        int index = 0, callStackPtr = 0;
        bool popping = false;

        if (!noLoad && !quiting && File.Exists(checkpointName)) {
            try {
                var lines = File.ReadAllLines(checkpointName);
                int lineNo = 0;
                var savedHeader = lines[lineNo++];
                if (savedHeader != filterSetHeader)
                    throw new InvalidDataException("Checkpoint filter set mismatch");
                var ss = lines[lineNo++].Split(',');
                for (int i = 0; i < ss.Length; i++) byteBoard[i] = byte.Parse(ss[i]);
                ss = lines[lineNo++].Split(',');
                for (int i = 0; i < ss.Length; i++) refStack[i] = int.Parse(ss[i]);
                ss = lines[lineNo++].Split(',');
                for (int i = 0; i < ss.Length; i++) callStack[i] = int.Parse(ss[i]);
                ss = lines[lineNo++].Split(',');
                depth = int.Parse(ss[0]);
                stackPtr = int.Parse(ss[1]);
                stackTopOriginal = int.Parse(ss[2]);
                stackLimit = int.Parse(ss[3]);
                index = int.Parse(ss[4]);
                popping = bool.Parse(ss[5]);
                callStackPtr = int.Parse(ss[6]);
                ss = lines[lineNo].Split(',');
                for (int i = 0; i < ss.Length; i++) filterCounts[i] = Num.Parse(ss[i]);
            } catch {
                Array.Clear(byteBoard);
                Array.Fill(byteBoard, (byte)255, Z + 1, boardSize - Z - 1);
                Array.Clear(refStack);
                refStack[0] = Z;
                Array.Clear(callStack);
                Array.Clear(filterCounts);
                depth = N; stackPtr = 1; stackTopOriginal = 1;
                stackLimit = stackSize; index = 0; callStackPtr = 0; popping = false;
            }
        }

        var swCheckpoint = System.Diagnostics.Stopwatch.StartNew();
        long lastSpinMs = 0;
        const string spinner = "⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏";
        int spinIdx = 0;
        ushort iterCounter = 0;

        for (;;) {
            if (++iterCounter == 0) {
                long elapsedMs = swCheckpoint.ElapsedMilliseconds;
                if (quiting || (!noSave && elapsedMs >= 60000L)) {
                    SaveCheckpoint(checkpointName, filterSetHeader, byteBoard, refStack, callStack,
                        depth, stackPtr, stackTopOriginal, stackLimit, index, popping,
                        callStackPtr, filterCounts);
                    swCheckpoint.Restart();
                    elapsedMs = 0;
                    lastSpinMs = 0;
                    if (quiting) return;
                }
                if (!quiet && elapsedMs - lastSpinMs >= 250) {
                    lastSpinMs = elapsedMs;
                    Console.Write($"\r{spinner[spinIdx++ % spinner.Length]}");
                }
            } else if (quiting) {
                SaveCheckpoint(checkpointName, filterSetHeader, byteBoard, refStack, callStack,
                    depth, stackPtr, stackTopOriginal, stackLimit, index, popping,
                    callStackPtr, filterCounts);
                return;
            }

            bool looping = stackPtr != 0;
            if (!popping && looping) {
                index = refStack[--stackPtr];

                int stackTopInner = stackPtr;
                if (++byteBoard[index - Z] == 0) refStack[stackTopInner++] = index - Z;
                if (++byteBoard[index - Y] == 0) refStack[stackTopInner++] = index - Y;
                if (++byteBoard[index - X] == 0) refStack[stackTopInner++] = index - X;
                if (++byteBoard[index + X] == 0) refStack[stackTopInner++] = index + X;
                if (++byteBoard[index + Y] == 0) refStack[stackTopInner++] = index + Y;
                if (++byteBoard[index + Z] == 0) refStack[stackTopInner++] = index + Z;

                if (depth == FilterDepth + 1) {
                    DispatchFiltersUnsafePinnedWork5(byteBoard, refStack, stackTopInner, stackLimit,
                        filterSet, filterCounts);
                } else if (depth > FilterDepth + 1) {
                    callStack[callStackPtr++] = stackPtr;
                    callStack[callStackPtr++] = stackLimit;
                    callStack[callStackPtr++] = stackTopOriginal;
                    callStack[callStackPtr++] = index;
                    depth--;
                    stackTopOriginal = stackPtr = stackTopInner;
                    continue;
                }
            }

            if (popping || looping) {
                --byteBoard[index - Z];
                --byteBoard[index - Y];
                --byteBoard[index - X];
                --byteBoard[index + X];
                --byteBoard[index + Y];
                --byteBoard[index + Z];
                refStack[--stackLimit] = index;
            }
            popping = false;

            if (stackPtr == 0) {
                Array.Copy(refStack, stackLimit, refStack, 0, stackTopOriginal);
                if (callStackPtr == 0) break;
                index = callStack[--callStackPtr];
                stackTopOriginal = callStack[--callStackPtr];
                stackLimit = callStack[--callStackPtr];
                stackPtr = callStack[--callStackPtr];
                depth++;
                popping = true;
            }
        }

        if (!quiting) {
            foreach (int f in filterSet)
                onFilterComplete(f, filterCounts[f]);
            try { if (File.Exists(checkpointName)) File.Delete(checkpointName); } catch { }
        }
    }
}










