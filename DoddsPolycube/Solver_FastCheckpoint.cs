namespace DoddsPolycube;

using Num = ulong;

// Fast checkpoint-resumable multi-filter computation
// ====================================================
//
// Same checkpointable outer loop as Exp6 (iterative, managed arrays, serializable),
// but when we reach FilterDepth+1, we pin the managed byteBoard and build a
// temporary byte** refStack to run DispatchFilters + CountExtensionsInner using
// fast unsafe pointer arithmetic — identical to Exp5's hot path.
//
// This gives us:
// - Checkpoint/resume for the outer loop (survives restarts)
// - Full unsafe pointer speed for the inner loop (no bounds checking)
// - The pin + pointer conversion at FilterDepth+1 is negligible overhead
//   (~thousands of times) vs the billions of inner loop iterations.

public partial class Program {

    private static Num[] CountExtensionsSubsetCheckpointedFast(int numPaths,
        bool quiet, bool noSave, bool skipLoad) {
        int totalFilters = MaxLeftStackLen + 1;
        Num[] counts = new Num[totalFilters];

        // Load previously completed filter results
        int loadedCount = 0;
        HashSet<int> alreadyDone = [];
        if (!skipLoad) {
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
                    } catch { /* corrupted file, recompute */ }
                }
            }
        }
        if (loadedCount > 0 && !quiet)
            Console.WriteLine($"  loaded {loadedCount} previously saved filter results");

        var remaining = new List<int>();
        for (int f = 0; f < totalFilters; f++)
            if (!alreadyDone.Contains(f))
                remaining.Add(f);

        if (remaining.Count == 0 || quiting) return counts;

        // Partition remaining filters into groups, interleaved for load balance
        int actualPaths = Math.Min(numPaths, remaining.Count);
        var groups = new HashSet<int>[actualPaths];
        for (int p = 0; p < actualPaths; p++)
            groups[p] = [];
        for (int i = 0; i < remaining.Count; i++)
            groups[i % actualPaths].Add(remaining[i]);

        // Progress tracking
        int completedFilters = loadedCount;
        var swProgress = System.Diagnostics.Stopwatch.StartNew();
        long lastProgressMs = 0;

        void OnFilterComplete(int filterIdx, Num count) {
            counts[filterIdx] = count;
            int done = Interlocked.Increment(ref completedFilters);

            if (!noSave) {
                var filename = $"trivial_{N}_{MaxLeftStackLen}_{filterIdx}.txt";
                var s = $"#{filterIdx} count={count:N0}";
                try {
                    File.WriteAllText(filename,
                        $"{count} 00:00:00{Environment.NewLine}{s}");
                } catch { /* best-effort */ }
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

        // Run each group in parallel — each with its own checkpointed outer loop
        Parallel.Invoke(groups.Select((filterSet, pathIdx) => (Action)(() => {
            if (quiting) return;
            CheckpointedFastWorker(filterSet, pathIdx, actualPaths,
                OnFilterComplete, noSave, skipLoad, quiet);
        })).ToArray());

        if (!quiet)
            Console.WriteLine($"\r  [{completedFilters}/{totalFilters}] elapsed={swProgress.Elapsed:hh\\:mm\\:ss}  ");

        return counts;
    }

    /// <summary>
    /// Checkpointable outer loop (managed arrays) + fast unsafe inner dispatch.
    /// </summary>
    private static void CheckpointedFastWorker(
        HashSet<int> filterSet, int pathIdx, int numPaths,
        Action<int, Num> onFilterComplete,
        bool noSave, bool skipLoad, bool quiet) {
        if (quiting) return;

        string checkpointName = $"checkpoint_{N}_fd{FilterDepth}_{numPaths}p_{pathIdx}.txt";
        string filterSetHeader = string.Join(',', filterSet.Order());
        int boardSize = (N + 2) * Z;
        int stackSize = (N - 2) * 4;

        // Managed arrays for the checkpointable outer loop
        byte[] byteBoard = new byte[boardSize];
        int[] refStack = new int[stackSize];
        int[] callStack = new int[stackSize];
        Num[] filterCounts = new Num[MaxLeftStackLen + 1];

        // Initialize board
        Array.Fill(byteBoard, (byte)255, Z + 1, boardSize - Z - 1);
        refStack[0] = Z;

        int depth = N, stackPtr = 1, stackTopOriginal = 1, stackLimit = stackSize;
        int index = 0, callStackPtr = 0;
        bool popping = false;

        // Try to load checkpoint
        if (!skipLoad && !quiting && File.Exists(checkpointName)) {
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

        // Main iterative outer loop (depths N → FilterDepth+1)
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
                    char spin = spinner[spinIdx++ % spinner.Length];
                    Console.Write($"\r{spin}");
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

                // Wind
                int stackTopInner = stackPtr;
                if (++byteBoard[index - Z] == 0) refStack[stackTopInner++] = index - Z;
                if (++byteBoard[index - Y] == 0) refStack[stackTopInner++] = index - Y;
                if (++byteBoard[index - X] == 0) refStack[stackTopInner++] = index - X;
                if (++byteBoard[index + X] == 0) refStack[stackTopInner++] = index + X;
                if (++byteBoard[index + Y] == 0) refStack[stackTopInner++] = index + Y;
                if (++byteBoard[index + Z] == 0) refStack[stackTopInner++] = index + Z;

                if (depth == FilterDepth + 1) {
                    // Pin the managed board and dispatch filters using fast unsafe pointers.
                    // The pin is done once per FilterDepth+1 node (~thousands of times),
                    // while the inner loop runs billions of times — negligible overhead.
                    DispatchFiltersUnsafePinned(byteBoard, refStack, stackTopInner, stackLimit,
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
                // Unwind
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

    /// <summary>
    /// Pins the managed byteBoard, converts int[] offsets to byte** pointers,
    /// and runs the fast unsafe DispatchFilters + CountExtensionsInner.
    /// </summary>
    private static unsafe void DispatchFiltersUnsafePinned(
        byte[] byteBoard, int[] refStackArr, int stackTopInnerIdx, int stackLimitIdx,
        HashSet<int> filterSet, Num[] filterCounts) {

        fixed (byte* boardBase = byteBoard) {
            // Build a temporary pointer-based stack on the native stack.
            // Only need the left side (0..stackTopInnerIdx) for the iteration.
            // Right side starts empty at the end for the inner recursion to use.
            int stackSize = (N - 2) * 4;
            byte** ptrStack = stackalloc byte*[stackSize];

            // Convert left-side int offsets to byte* pointers
            for (int i = 0; i < stackTopInnerIdx; i++)
                ptrStack[i] = boardBase + refStackArr[i];

            byte** ptrBase = ptrStack;
            byte** stackTop1 = ptrStack + stackTopInnerIdx;
            byte** stackTop2 = ptrStack + stackLimitIdx; // right stack starts at limit

            // === FilterDepth iteration: same as Exp5's DispatchFilters ===
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

                    filterCounts[filterIdx] += CountExtInner(
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
            // Restore: copy right back to left in the pointer stack
            while (stackTop1 != stackTopOriginal)
                *stackTop1++ = *stackTop2++;

            // Copy the restored left side back to the managed int[] array,
            // so the outer loop's managed state stays consistent.
            for (int i = 0; i < stackTopInnerIdx; i++)
                refStackArr[i] = (int)(ptrStack[i] - boardBase);
        }
    }

    /// <summary>
    /// Fast unsafe inner recursion — identical to Exp5's CountExtensionsInner.
    /// Standalone static method (not a local function) so it can be called from Exp8.
    /// </summary>
    private static unsafe Num CountExtInner(
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

            if (depth == 4) {
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
            } else {
                count += CountExtInner(refStackBase, depth - 1, stackTopInner, stackTop2);
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
}




