namespace DoddsPolycube;

using Num = ulong;

// Checkpoint-resumable multi-filter single-pass optimization
// ============================================================
//
// Combines the Exp5 approach (partition filters into groups, traverse shared tree once per
// group) with the CountExtensionsSubsetStack approach (iterative outer loop with explicit
// call stack that can be serialized/deserialized for checkpoint/resume).
//
// Architecture:
// - OUTER loop (depths N → FilterDepth+1): iterative with explicit call stack.
//   Checkpoints are written here (every ~30s or on Ctrl+C). State is small (~10KB).
// - FILTER level (depth FilterDepth): iterates stack positions, dispatches to inner loop.
// - INNER loop (depths FilterDepth-1 → 4 + Work4): recursive + unsafe pointers for
//   maximum performance. Runs to completion without checkpoint checks.
//
// Checkpoint file format (text, one line per section):
//   Line 0: byteBoard as comma-separated bytes
//   Line 1: refStack as comma-separated ints
//   Line 2: callStack as comma-separated ints
//   Line 3: scalar state (depth,stackPtr,stackTopOriginal,stackLimit,index,popping,callStackPtr)
//   Line 4: per-filter counts as comma-separated ulongs
//
// For N=23: board is (23+2)*Z bytes, refStack is (23-2)*4 ints, callStack is (23-2)*4 ints.
// Total checkpoint size: ~20KB text.

public partial class Program {
    /// <summary>
    /// Like CountExtensionsSubsetAllFiltersParallel, but each path uses an iterative outer
    /// loop with periodic checkpointing so progress survives process restarts.
    /// </summary>
    private static Num[] CountExtensionsSubsetCheckpointed(int numPaths,
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
        pathCounts = new Num[actualPaths];
        Parallel.Invoke(groups.Select((filterSet, pathIdx) => (Action)(() => {
            if (quiting) return;
            CountExtensionsSubsetFiltersCheckpointed(filterSet, pathIdx, actualPaths,
                OnFilterComplete, noSave, skipLoad, quiet);
        })).ToArray());

        if (!quiet)
            Console.WriteLine($"\r  [{completedFilters}/{totalFilters}] elapsed={swProgress.Elapsed:hh\\:mm\\:ss}  ");

        return counts;
    }

    private static Num[] pathCounts;

    /// <summary>
    /// Processes a SET of filter values using an iterative outer loop (depths N → FilterDepth+1)
    /// with periodic checkpointing. The inner loop (depths FilterDepth-1 → 4 + Work4) remains
    /// recursive+unsafe for performance.
    /// </summary>
    private static void CountExtensionsSubsetFiltersCheckpointed(
        HashSet<int> filterSet, int pathIdx, int numPaths,
        Action<int, Num> onFilterComplete,
        bool noSave, bool skipLoad, bool quiet) {
        if (quiting) return;

        string checkpointName = $"checkpoint_{N}_fd{FilterDepth}_{numPaths}p_{pathIdx}.txt";
        // Sorted filter set string used as header for validation
        string filterSetHeader = string.Join(',', filterSet.Order());
        int boardSize = (N + 2) * Z;
        int stackSize = (N - 2) * 4;

        // Use managed arrays for the outer iterative loop (serializable)
        byte[] byteBoard = new byte[boardSize];
        int[] refStack = new int[stackSize];
        int[] callStack = new int[stackSize]; // 4 ints per frame, max depth frames
        Num[] filterCounts = new Num[MaxLeftStackLen + 1];

        // Initialize board: first Z+1 bytes are 0, rest 255
        Array.Fill(byteBoard, (byte)255, Z + 1, boardSize - Z - 1);
        refStack[0] = Z; // seed with first board position

        // Iterative state
        int depth = N, stackPtr = 1, stackTopOriginal = 1, stackLimit = stackSize;
        int index = 0, callStackPtr = 0;
        bool popping = false;

        // Try to load checkpoint
        if (!skipLoad && !quiting && File.Exists(checkpointName)) {
            try {
                var lines = File.ReadAllLines(checkpointName);
                int lineNo = 0;
                // Line 0: filter set header — must match exactly
                var savedHeader = lines[lineNo++];
                if (savedHeader != filterSetHeader)
                    throw new InvalidDataException(
                        $"Checkpoint filter set mismatch: expected [{filterSetHeader}], got [{savedHeader}]");
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
                // Corrupted checkpoint — start fresh
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
            // Only check timers every 65536 iterations to avoid Stopwatch overhead
            if (++iterCounter == 0) {
                long elapsedMs = swCheckpoint.ElapsedMilliseconds;

                // Checkpoint save: on quit signal or periodic timer (~1m)
                if (quiting || (!noSave && elapsedMs >= 60000L)) {
                    SaveCheckpoint(checkpointName, filterSetHeader, byteBoard, refStack, callStack,
                        depth, stackPtr, stackTopOriginal, stackLimit, index, popping,
                        callStackPtr, filterCounts);
                    swCheckpoint.Restart();
                    elapsedMs = 0;
                    lastSpinMs = 0;
                    if (quiting) return;
                    Num totalCounts = 0;
                    foreach (Num fc in filterCounts) totalCounts += fc;
                    pathCounts[pathIdx] = totalCounts;
                    Num grandTotalCounts = 0;
                    foreach (Num pc in pathCounts) grandTotalCounts += pc;
                    if (!quiet)
                        Console.WriteLine($"\rpath[{pathIdx}]={totalCounts:N0}, total={grandTotalCounts:N0}");
                }

                // Spinner update (~250ms)
                if (!quiet && elapsedMs - lastSpinMs >= 250) {
                    lastSpinMs = elapsedMs;
                    char spin = spinner[spinIdx++ & 7];
                    Console.Write($"\r{spin}");
                }
            } else if (quiting) {
                // Always check quit flag immediately
                SaveCheckpoint(checkpointName, filterSetHeader, byteBoard, refStack, callStack,
                    depth, stackPtr, stackTopOriginal, stackLimit, index, popping,
                    callStackPtr, filterCounts);
                return;
            }

            bool looping = stackPtr != 0, skipUnwind = false;
            if (!popping && looping) {
                index = refStack[--stackPtr];

                // Wind: increment 6 neighbors, push zero-valued ones
                int stackTopInner = stackPtr;
                if (++byteBoard[index - Z] == 0) refStack[stackTopInner++] = index - Z;
                if (++byteBoard[index - Y] == 0) refStack[stackTopInner++] = index - Y;
                if (++byteBoard[index - X] == 0) refStack[stackTopInner++] = index - X;
                if (++byteBoard[index + X] == 0) refStack[stackTopInner++] = index + X;
                if (++byteBoard[index + Y] == 0) refStack[stackTopInner++] = index + Y;
                if (++byteBoard[index + Z] == 0) refStack[stackTopInner++] = index + Z;

                if (depth == FilterDepth + 1) {
                    // We're one level above FilterDepth. Dispatch filters using the
                    // fast unsafe inner path. This runs to completion (no checkpointing
                    // inside), but each subtree takes at most seconds.
                    DispatchFiltersSafe(byteBoard, refStack, stackTopInner, stackLimit,
                        filterSet, filterCounts);
                } else if (depth > FilterDepth + 1) {
                    // Push frame and recurse (iteratively)
                    callStack[callStackPtr++] = stackPtr;
                    callStack[callStackPtr++] = stackLimit;
                    callStack[callStackPtr++] = stackTopOriginal;
                    callStack[callStackPtr++] = index;
                    depth--;
                    stackTopOriginal = stackPtr = stackTopInner;
                    continue; // skip unwind — we're "entering" the recursive call
                }
                // If depth == FilterDepth + 1, we fall through to unwind after dispatch
            }

            if (popping || looping) {
                if (!skipUnwind) {
                    // Unwind: decrement 6 neighbors
                    --byteBoard[index - Z];
                    --byteBoard[index - Y];
                    --byteBoard[index - X];
                    --byteBoard[index + X];
                    --byteBoard[index + Y];
                    --byteBoard[index + Z];
                }
                refStack[--stackLimit] = index;
            }
            popping = false;

            if (stackPtr == 0) {
                // Restore stack: copy right side back to left
                Array.Copy(refStack, stackLimit, refStack, 0, stackTopOriginal);

                if (callStackPtr == 0) break; // done!
                // Pop frame
                index = callStack[--callStackPtr];
                stackTopOriginal = callStack[--callStackPtr];
                stackLimit = callStack[--callStackPtr];
                stackPtr = callStack[--callStackPtr];
                depth++;
                popping = true; // next iteration does unwind
            }
        }

        // Completed — report results and clean up checkpoint
        if (!quiting) {
            foreach (int f in filterSet)
                onFilterComplete(f, filterCounts[f]);
            try { if (File.Exists(checkpointName)) File.Delete(checkpointName); } catch { }
        }
    }

    /// <summary>
    /// At depth FilterDepth: iterate stack positions, dispatch each matching filter's subtree
    /// using the index-based inner recursion on managed arrays.
    /// </summary>
    private static void DispatchFiltersSafe(byte[] byteBoard, int[] refStackArr,
        int stackTopInnerIdx, int stackLimitIdx,
        HashSet<int> filterSet, Num[] filterCounts) {
        int stackPtr = stackTopInnerIdx;
        int stackLimit = stackLimitIdx;
        int stackTopOriginal = stackPtr;

        while (stackPtr != 0) {
            int indexOffset = refStackArr[--stackPtr];
            int filterIdx = stackPtr;
            if (filterSet.Contains(filterIdx)) {
                int stackTopInner = stackPtr;

                if (++byteBoard[indexOffset - Z] == 0) refStackArr[stackTopInner++] = indexOffset - Z;
                if (++byteBoard[indexOffset - Y] == 0) refStackArr[stackTopInner++] = indexOffset - Y;
                if (++byteBoard[indexOffset - X] == 0) refStackArr[stackTopInner++] = indexOffset - X;
                if (++byteBoard[indexOffset + X] == 0) refStackArr[stackTopInner++] = indexOffset + X;
                if (++byteBoard[indexOffset + Y] == 0) refStackArr[stackTopInner++] = indexOffset + Y;
                if (++byteBoard[indexOffset + Z] == 0) refStackArr[stackTopInner++] = indexOffset + Z;

                filterCounts[filterIdx] += CountExtensionsInnerSafe(byteBoard, refStackArr,
                    FilterDepth - 1, stackTopInner, stackLimit);

                --byteBoard[indexOffset - Z];
                --byteBoard[indexOffset - Y];
                --byteBoard[indexOffset - X];
                --byteBoard[indexOffset + X];
                --byteBoard[indexOffset + Y];
                --byteBoard[indexOffset + Z];
            }

            refStackArr[--stackLimit] = indexOffset;
        }

        while (stackPtr != stackTopOriginal)
            refStackArr[stackPtr++] = refStackArr[stackLimit++];
    }

    /// <summary>
    /// Index-based inner recursion (managed arrays). For the hot path below FilterDepth.
    /// Uses the same algorithm as CountExtensionsInner but with int[] indices instead of byte**.
    /// </summary>
    private static Num CountExtensionsInnerSafe(byte[] board, int[] stack,
        int depth, int stackPtr, int stackLimit) {
        Num count = 0;
        int stackTopOriginal = stackPtr;
        while (stackPtr != 0) {
            int index = stack[--stackPtr];
            int stackTopInner = stackPtr;

            if (++board[index - Z] == 0) stack[stackTopInner++] = index - Z;
            if (++board[index - Y] == 0) stack[stackTopInner++] = index - Y;
            if (++board[index - X] == 0) stack[stackTopInner++] = index - X;
            if (++board[index + X] == 0) stack[stackTopInner++] = index + X;
            if (++board[index + Y] == 0) stack[stackTopInner++] = index + Y;
            if (++board[index + Z] == 0) stack[stackTopInner++] = index + Z;

            if (depth == 4) {
                int stackTop = stackTopInner, stackTopTemp = stackTopInner;
                int length = stackTopInner, lengthPlus = (stackTopInner << 1) - 511;
                count += (Num)(length * (length - 1) * (length - 2) / 6);
                for (; stackTopTemp != 0;) {
                    int i = stack[--stackTopTemp], neighbours = 0, subCount = 128, localCount = 0;
                    byte v;
                    int ii;
                    if ((v = board[ii = i - Z]) > 127) {
                        localCount += board[ii] = --v;
                        subCount += board[ii - Z] + board[ii - X] + board[ii - Y] + board[ii + X] + board[ii + Y];
                        neighbours++;
                    }
                    if ((v = board[ii = i - Y]) > 127) {
                        localCount += board[ii] = --v;
                        subCount += board[ii - Y] + board[ii - X] + board[ii - Z] + board[ii + X] + board[ii + Z];
                        neighbours++;
                    }
                    if ((v = board[ii = i - X]) > 127) {
                        localCount += board[ii] = --v;
                        subCount += board[ii - X] + board[ii - Y] + board[ii - Z] + board[ii + Y] + board[ii + Z];
                        neighbours++;
                    }
                    if ((v = board[ii = i + X]) > 127) {
                        localCount += board[ii] = --v;
                        subCount += board[ii + X] + board[ii + Y] + board[ii + Z] + board[ii - Y] + board[ii - Z];
                        neighbours++;
                    }
                    if ((v = board[ii = i + Y]) > 127) {
                        localCount += board[ii] = --v;
                        subCount += board[ii + Y] + board[ii + X] + board[ii + Z] + board[ii - X] + board[ii - Z];
                        neighbours++;
                    }
                    if ((v = board[ii = i + Z]) > 127) {
                        localCount += board[ii] = --v;
                        subCount += board[ii + Z] + board[ii + X] + board[ii + Y] + board[ii - X] + board[ii - Y];
                        neighbours++;
                    }
                    count += (Num)localCount +
                             (Num)((subCount >> 8) + (neighbours * (neighbours + lengthPlus) >> 1));
                }
                while (stackTop != 0) {
                    int i = stack[--stackTop];
                    int ii;
                    board[ii = i - Z] |= (byte)(board[ii] >> 4);
                    board[ii = i - Y] |= (byte)(board[ii] >> 4);
                    board[ii = i - X] |= (byte)(board[ii] >> 4);
                    board[ii = i + X] |= (byte)(board[ii] >> 4);
                    board[ii = i + Y] |= (byte)(board[ii] >> 4);
                    board[ii = i + Z] |= (byte)(board[ii] >> 4);
                }
            } else {
                count += CountExtensionsInnerSafe(board, stack, depth - 1, stackTopInner, stackLimit);
            }

            --board[index - Z];
            --board[index - Y];
            --board[index - X];
            --board[index + X];
            --board[index + Y];
            --board[index + Z];

            stack[--stackLimit] = index;
        }
        while (stackPtr != stackTopOriginal)
            stack[stackPtr++] = stack[stackLimit++];
        return count;
    }

    private static void SaveCheckpoint(string filename, string filterSetHeader,
        byte[] byteBoard, int[] refStack,
        int[] callStack, int depth, int stackPtr, int stackTopOriginal, int stackLimit,
        int index, bool popping, int callStackPtr, Num[] filterCounts) {
        try {
            var tmpName = filename + ".tmp";
            using (var w = new StreamWriter(tmpName)) {
                w.WriteLine(filterSetHeader);
                w.WriteLine(string.Join(',', byteBoard));
                w.WriteLine(string.Join(',', refStack));
                w.WriteLine(string.Join(',', callStack));
                w.WriteLine($"{depth},{stackPtr},{stackTopOriginal},{stackLimit},{index},{popping},{callStackPtr}");
                w.WriteLine(string.Join(',', filterCounts));
            }
            // Atomic rename to avoid corrupted checkpoints
            if (File.Exists(filename)) File.Delete(filename);
            File.Move(tmpName, filename);
        } catch {
            // best-effort checkpoint
        }

        if (Program.checkpointLog) {
            try {
                var logName = $"checkpoint_{N}_fd{FilterDepth}.log";
                Num runningTotal = 0;
                int completedFilters = 0;
                for (int i = 0; i < filterCounts.Length; i++) {
                    runningTotal += filterCounts[i];
                    if (filterCounts[i] > 0) completedFilters++;
                }
                var line = $"{DateTime.UtcNow:O} {filename}" +
                           $" depth={depth} stackPtr={stackPtr} callStackPtr={callStackPtr}" +
                           $" filters={completedFilters}/{filterCounts.Length}" +
                           $" runningTotal={runningTotal}" +
                           (quiting ? " QUITING" : "");
                File.AppendAllText(logName, line + Environment.NewLine);
            } catch {
                // best-effort logging
            }
        }
    }
}










