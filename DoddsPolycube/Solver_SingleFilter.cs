namespace DoddsPolycube;

using Num = ulong;

public partial class Program {
    private static unsafe Num CountExtensionsSubsetUnsafe2(int filter) {
        if (quiting) return 0;
        // could use ints or shorts as offsets to save memory, but it's faster to directly store the
        // pointers to avoid adding pointer offsets at every lookup
        byte* byteBoard = stackalloc byte[(N + 2) * Z];
        // total length of the two stacks is at most 4n-9. One stack grows from the left, the other
        // stack grows from the right
        byte** refStack = stackalloc byte*[(N - 2) * 4];
        // seeded with first index of the byte board as the only allowed extension
        *refStack = byteBoard += Z;

        for (byte* i = byteBoard + (N + 1) * Z; --i != byteBoard;)
            // the first Z + 1 bytes are disallowed extensions; first Z are less than the minimum,
            // last 1 due to edge case of initial polycube having no neighbours
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
                        // if multithreading is not wanted, remove "if (condition)" from this else statement
                        count += CountExtensions(depth - 1, stackTopInner, stackTop2);
                    }

                    --*(index - Z);
                    --*(index - Y);
                    --*(index - X);
                    --*(index + X);
                    --*(index + Y);
                    --*(index + Z);
                }

                // doing this push before the recursion would add one extra unnecessary element to the stack at
                // each level of recursion
                *--stackTop2 = index;
            }
            while (stackTop1 != stackTopOriginal)
                *stackTop1++ = *stackTop2++;
            return count;
        }
    }
}