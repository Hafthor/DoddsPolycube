namespace DoddsPolycube;

using Num = ulong;

// Two-phase split: single-threaded precalculation + parallel execution
// =====================================================================
//
// Phase A (single-threaded):
//   Traverse the top of the recursion tree (depths N → splitDepth) once,
//   collecting "work items" — board+stack snapshots at splitDepth nodes.
//   This is very fast since the tree above splitDepth is tiny.
//
// Phase B (parallel, fast unsafe):
//   Process all work items in parallel via Parallel.ForEach. Each work item
//   gets its own unsafe byte* board and byte** refStack, recurses down to
//   FilterDepth, dispatches all filters, and runs CountExtensionsInner at
//   full speed with no bounds checking.
//
// Note: The shared tree above splitDepth is negligible compared to the work
//   below it. The main advantage of this approach over Exp5 is:
//   - Better load balancing: thousands of fine-grained work items vs P coarse
//     paths, letting Parallel.ForEach distribute work more evenly.
//   - Combines unsafe pointer speed (Exp5) with the ability to save/load
//     per-filter results (like Exp6), without the managed-array overhead.
//   - Simpler parallelism: no need to partition filter sets across paths.

public partial class Program {

    // A work item captures the board + stack state at a FilterDepth+1 node.
    // The board is shared but mutated, so we must snapshot it.
    private struct FilterWorkItem {
        public byte[] Board;       // full board snapshot
        public int[] LeftStack;    // refStack[0..StackTop) — the active positions
        public int StackTop;       // number of elements in left stack
        public int Depth;          // depth at which this item was captured
    }

    /// <summary>
    /// Two-phase approach: single-threaded shared tree + parallel filter dispatch.
    /// </summary>
    private static Num CountExtensionsTwoPhase(int filter) {
        if (filter != 0) return 0;
        int totalFilters = MaxLeftStackLen + 1;
        Num[] counts = new Num[totalFilters];

        // Load previously completed filter results
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
                    } catch { /* corrupted file, recompute */ }
                }
            }
        }
        if (loadedCount > 0 && !quiet)
            Console.WriteLine($"  loaded {loadedCount} previously saved filter results");

        Num total = 0;
        if (alreadyDone.Count == totalFilters || quiting) {
            foreach (Num count in counts) total += count;
            return total;
        }

        // ---- Phase A: single-threaded shared tree traversal ----
        // Choose split depth: lower = more work items (finer parallelism).
        // Target ~1000s of items. Empirically, depth N-5 gives ~1000-10000 items.
        // Must be >= FilterDepth+1 so each item can dispatch filters.
        int splitDepth = Math.Max(FilterDepth + 1, N - 5);
        if (!quiet)
            Console.WriteLine($"  Phase A: traversing shared tree to depth {splitDepth} (single-threaded)...");

        var swA = System.Diagnostics.Stopwatch.StartNew();
        var workItems = CollectWorkItems(splitDepth, quiet);
        swA.Stop();

        if (quiting) {
            foreach (Num count in counts) total += count;
            return total;
        }
        if (!quiet)
            Console.WriteLine($"\r  Phase A complete: {workItems.Count:N0} work items in {swA.Elapsed:hh\\:mm\\:ss}    ");

        // ---- Phase B: parallel filter dispatch ----
        if (!quiet)
            Console.WriteLine($"  Phase B: processing {workItems.Count:N0} work items across {numPaths} threads...");

        var swB = System.Diagnostics.Stopwatch.StartNew();
        int completedItems = 0;
        long lastProgressMs = 0;

        // Each work item processes ALL non-done filters, accumulating into thread-local counts
        // then merged. Using Parallel.ForEach with thread-local accumulators.
        var lockObj = new object();
        Parallel.ForEach(
            workItems,
            new ParallelOptions { MaxDegreeOfParallelism = numPaths },
            () => new Num[totalFilters], // thread-local accumulator
            (workItem, state, localCounts) => {
                if (quiting) { state.Stop(); return localCounts; }
                ProcessWorkItemUnsafe(workItem, alreadyDone, localCounts);

                // Progress
                int done = Interlocked.Increment(ref completedItems);
                if (!quiet) {
                    long nowMs = swB.ElapsedMilliseconds;
                    long last = Interlocked.Read(ref lastProgressMs);
                    if (nowMs - last >= 1000 &&
                        Interlocked.CompareExchange(ref lastProgressMs, nowMs, last) == last) {
                        Console.Write($"\r  [{done:N0}/{workItems.Count:N0}] elapsed={swB.Elapsed:hh\\:mm\\:ss}  ");
                    }
                }
                return localCounts;
            },
            localCounts => {
                lock (lockObj) {
                    for (int f = 0; f < totalFilters; f++)
                        counts[f] += localCounts[f];
                }
            }
        );
        swB.Stop();

        if (!quiet)
            Console.WriteLine($"\r  Phase B complete: {swB.Elapsed:hh\\:mm\\:ss}                          ");

        // Save per-filter results
        if (!noSave && !quiting) {
            for (int f = 0; f < totalFilters; f++) {
                if (alreadyDone.Contains(f)) continue;
                var filename = $"trivial_{N}_{MaxLeftStackLen}_{f}.txt";
                var s = $"#{f} count={counts[f]:N0}";
                try {
                    File.WriteAllText(filename,
                        $"{counts[f]} 00:00:00{Environment.NewLine}{s}");
                } catch { /* best-effort */ }
            }
        }

        total = 0;
        foreach (Num count in counts) total += count;
        return total;
    }

    /// <summary>
    /// Phase A: traverse the shared tree (depths N → splitDepth) exactly once,
    /// collecting a work item at each splitDepth node.
    /// splitDepth must be >= FilterDepth+1 (so each work item can dispatch filters).
    /// Lower splitDepth = more work items (finer parallelism, more memory).
    /// Higher splitDepth = fewer work items (coarser parallelism, less memory).
    /// </summary>
    private static List<FilterWorkItem> CollectWorkItems(int splitDepth, bool quiet) {
        int boardSize = (N + 2) * Z;
        int stackSize = (N - 2) * 4;

        byte[] byteBoard = new byte[boardSize];
        int[] refStack = new int[stackSize];
        int[] callStack = new int[stackSize];

        Array.Fill(byteBoard, (byte)255, Z + 1, boardSize - Z - 1);
        refStack[0] = Z;

        int depth = N, stackPtr = 1, stackTopOriginal = 1, stackLimit = stackSize;
        int index = 0, callStackPtr = 0;
        bool popping = false;

        var workItems = new List<FilterWorkItem>();

        // Try to load checkpoint (includes partially collected work items)
        // For simplicity, we don't checkpoint the work item list — Phase A is fast
        // (single traversal of shared tree). For N=23 this takes minutes, not hours.
        // If interrupted, it restarts Phase A from scratch.

        var sw = System.Diagnostics.Stopwatch.StartNew();
        long lastSpinMs = 0;
        const string spinner = "⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏";
        int spinIdx = 0;
        ushort iterCounter = 0;

        for (;;) {
            if (++iterCounter == 0) {
                long elapsedMs = sw.ElapsedMilliseconds;

                if (quiting) return workItems;

                if (!quiet && elapsedMs - lastSpinMs >= 250) {
                    lastSpinMs = elapsedMs;
                    char spin = spinner[spinIdx++ % spinner.Length];
                    Console.Write($"\r  {spin} collecting work items: {workItems.Count:N0} (depth={depth})  ");
                }
            } else if (quiting) {
                return workItems;
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

                if (depth == splitDepth) {
                    // Capture work item: snapshot board and left stack
                    var item = new FilterWorkItem {
                        Board = (byte[])byteBoard.Clone(),
                        LeftStack = new int[stackTopInner],
                        StackTop = stackTopInner,
                        Depth = depth,
                    };
                    Array.Copy(refStack, 0, item.LeftStack, 0, stackTopInner);
                    workItems.Add(item);
                    // Fall through to unwind
                } else if (depth > splitDepth) {
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

        return workItems;
    }

    /// <summary>
    /// Phase B: process a single work item using fast unsafe pointers.
    /// If the work item was captured above FilterDepth+1, recurses down through
    /// the remaining shared tree levels before dispatching filters.
    /// </summary>
    private static unsafe void ProcessWorkItemUnsafe(FilterWorkItem item,
        HashSet<int> alreadyDone, Num[] filterCounts) {
        int boardSize = (N + 2) * Z;
        int stackSize = (N - 2) * 4;

        byte* byteBoard = stackalloc byte[boardSize];
        byte** refStack = stackalloc byte*[stackSize];

        for (int i = 0; i < boardSize; i++)
            byteBoard[i] = item.Board[i];

        // Restore left stack only — right stack starts empty at the end
        int stackTop = item.StackTop;
        for (int i = 0; i < stackTop; i++)
            refStack[i] = byteBoard + item.LeftStack[i];

        byte** refStackBase = refStack;
        byte** stackTop1 = refStack + stackTop;
        byte** stackTop2 = refStack + stackSize; // right stack starts empty

        Num[] localCounts = new Num[MaxLeftStackLen + 1];

        // The captured state is AFTER Wind at splitDepth. The left stack contains
        // the children ready for the next level (splitDepth - 1).
        // So we process as if we're already inside the loop at splitDepth,
        // having just wound, and need to recurse into depth splitDepth-1.
        if (item.Depth - 1 == FilterDepth) {
            // Next level IS FilterDepth — dispatch filters directly
            DispatchFiltersUnsafe(refStackBase, stackTop1, stackTop2, alreadyDone, localCounts);
        } else {
            // Recurse from depth-1 down to FilterDepth+1
            ProcessSubtreeUnsafe(refStackBase, item.Depth - 1, stackTop1, stackTop2, alreadyDone, localCounts);
        }

        // Merge into caller's filterCounts
        for (int f = 0; f <= MaxLeftStackLen; f++)
            filterCounts[f] += localCounts[f];
    }

    /// <summary>
    /// Recursive shared tree traversal from captured depth down to FilterDepth+1,
    /// then dispatches filters at FilterDepth.
    /// </summary>
    private static unsafe void ProcessSubtreeUnsafe(byte** refStackBase,
        int depth, byte** stackTop1, byte** stackTop2,
        HashSet<int> alreadyDone, Num[] filterCounts) {
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

            if (depth == FilterDepth + 1) {
                // Dispatch filters
                DispatchFiltersUnsafe(refStackBase, stackTopInner, stackTop2,
                    alreadyDone, filterCounts);
            } else {
                ProcessSubtreeUnsafe(refStackBase, depth - 1, stackTopInner, stackTop2,
                    alreadyDone, filterCounts);
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
    }

    /// <summary>
    /// Dispatch filters at FilterDepth, calling fast inner recursion for each.
    /// </summary>
    private static unsafe void DispatchFiltersUnsafe(byte** refStackBase,
        byte** stackTop1, byte** stackTop2,
        HashSet<int> alreadyDone, Num[] filterCounts) {
        byte** stackTopOriginal = stackTop1;
        while (stackTop1 != refStackBase) {
            byte* idx = *--stackTop1;
            int filterIdx = (int)(stackTop1 - refStackBase);
            if (filterIdx <= MaxLeftStackLen && !alreadyDone.Contains(filterIdx)) {
                byte** stackTopInner = stackTop1;

                if (++*(idx - Z) == 0) *stackTopInner++ = idx - Z;
                if (++*(idx - Y) == 0) *stackTopInner++ = idx - Y;
                if (++*(idx - X) == 0) *stackTopInner++ = idx - X;
                if (++*(idx + X) == 0) *stackTopInner++ = idx + X;
                if (++*(idx + Y) == 0) *stackTopInner++ = idx + Y;
                if (++*(idx + Z) == 0) *stackTopInner++ = idx + Z;

                filterCounts[filterIdx] += CountExtensionsInnerUnsafe(
                    refStackBase, FilterDepth - 1, stackTopInner, stackTop2);

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
    }

    /// <summary>
    /// Fast unsafe inner recursion (same as Exp5 CountExtensionsInner).
    /// </summary>
    private static unsafe Num CountExtensionsInnerUnsafe(
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
                count += CountExtensionsInnerUnsafe(refStackBase, depth - 1, stackTopInner, stackTop2);
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



















