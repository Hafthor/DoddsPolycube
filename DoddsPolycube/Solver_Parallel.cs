namespace DoddsPolycube;

using Num = ulong;

// Multi-filter single-pass optimization with parallel partitioning
// =================================================================
//
// Key insight: every filter value 0..MaxLeftStackLen traverses the IDENTICAL recursion tree
// at depths N down to FilterDepth+1. Instead of running all 43 filters as separate tasks
// (each redundantly traversing the shared prefix), we partition filters into P groups of
// ~43/P filters each. Each group traverses the shared prefix once and processes only its
// assigned filters at the FilterDepth gate. Groups run in parallel via Parallel.Invoke.
//
// Features:
// - Progress reporting: shared atomic counter tracks completed filter-subtrees
// - Save/Resume: completed per-filter results are saved to trivial_N_M_F.txt files;
//   on restart, already-saved filters are loaded and excluded from computation.

public partial class Program {
    /// <summary>
    /// Partitions filter values into numPaths groups and runs each group in parallel.
    /// Loads previously saved filter results and skips them. Saves results per-filter
    /// as each path completes. Reports progress via shared counter.
    /// </summary>
    private static Num[] CountExtensionsSubsetAllFiltersParallel(int numPaths,
        bool quiet, bool noSave, bool skipLoad) {
        int totalFilters = MaxLeftStackLen + 1; // 0..42 = 43 filters
        Num[] counts = new Num[totalFilters];

        // Load previously saved results
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
                    } catch {
                        // corrupted file, recompute
                    }
                }
            }
        }
        if (loadedCount > 0 && !quiet)
            Console.WriteLine($"  loaded {loadedCount} previously saved filter results");

        // Determine which filters still need computing
        var remaining = new List<int>();
        for (int f = 0; f < totalFilters; f++)
            if (!alreadyDone.Contains(f))
                remaining.Add(f);

        if (remaining.Count == 0) return counts;
        if (quiting) return counts;

        // Partition remaining filters into numPaths groups, interleaved for load balance
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

        // Callback invoked when a single filter subtree completes (from any path/thread)
        void OnFilterComplete(int filterIdx, Num count) {
            counts[filterIdx] = count;
            int done = Interlocked.Increment(ref completedFilters);

            if (!noSave) {
                var filename = $"trivial_{N}_{MaxLeftStackLen}_{filterIdx}.txt";
                var s = $"#{filterIdx} count={count:N0}";
                try {
                    File.WriteAllText(filename,
                        $"{count} 00:00:00{Environment.NewLine}{s}");
                } catch {
                    // best-effort save
                }
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

        // Run each group in parallel
        Parallel.Invoke(groups.Select(filterSet => (Action)(() => {
            if (quiting) return;
            CountExtensionsSubsetFilters(filterSet, OnFilterComplete);
        })).ToArray());

        if (!quiet)
            Console.WriteLine($"\r  [{completedFilters}/{totalFilters}] elapsed={swProgress.Elapsed:hh\\:mm\\:ss}  ");

        return counts;
    }

    /// <summary>
    /// Processes a specific SET of filter values in a single recursion pass.
    /// The tree above FilterDepth is traversed once; at FilterDepth only the
    /// specified filter positions are processed.
    /// Calls onFilterComplete(filterIdx, count) as each filter finishes.
    /// </summary>
    private static unsafe void CountExtensionsSubsetFilters(
        HashSet<int> filterSet, Action<int, Num> onFilterComplete) {
        if (quiting) return;
        byte* byteBoard = stackalloc byte[(N + 2) * Z];
        byte** refStack = stackalloc byte*[(N - 2) * 4];
        *refStack = byteBoard += Z;

        for (byte* i = byteBoard + (N + 1) * Z; --i != byteBoard;)
            *i = 255;

        Num[] counts = new Num[MaxLeftStackLen + 1];
        CountExtensionsOuter(N, refStack + 1, refStack + (N - 2) * 4);
        // Report final results for each filter in this set
        if (!quiting)
            foreach (int f in filterSet)
                onFilterComplete(f, counts[f]);

        // Depths N..FilterDepth+1: shared recursion tree, traversed once.
        void CountExtensionsOuter(int depth, byte** stackTop1, byte** stackTop2) {
            byte** stackTopOriginal = stackTop1;
            while (stackTop1 != refStack) {
                if (quiting) break;
                byte* index = *--stackTop1;
                byte** stackTopInner = stackTop1;

                if (++*(index - Z) == 0) *stackTopInner++ = index - Z;
                if (++*(index - Y) == 0) *stackTopInner++ = index - Y;
                if (++*(index - X) == 0) *stackTopInner++ = index - X;
                if (++*(index + X) == 0) *stackTopInner++ = index + X;
                if (++*(index + Y) == 0) *stackTopInner++ = index + Y;
                if (++*(index + Z) == 0) *stackTopInner++ = index + Z;

                if (depth == FilterDepth + 1) {
                    DispatchFilters(stackTopInner, stackTop2);
                } else {
                    CountExtensionsOuter(depth - 1, stackTopInner, stackTop2);
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

        // At FilterDepth: iterate stack positions, only process those in filterSet.
        void DispatchFilters(byte** stackTop1, byte** stackTop2) {
            byte** stackTopOriginal = stackTop1;
            while (stackTop1 != refStack) {
                byte* index = *--stackTop1;
                int filterIdx = (int)(stackTop1 - refStack);
                if (filterSet.Contains(filterIdx)) {
                    byte** stackTopInner = stackTop1;

                    if (++*(index - Z) == 0) *stackTopInner++ = index - Z;
                    if (++*(index - Y) == 0) *stackTopInner++ = index - Y;
                    if (++*(index - X) == 0) *stackTopInner++ = index - X;
                    if (++*(index + X) == 0) *stackTopInner++ = index + X;
                    if (++*(index + Y) == 0) *stackTopInner++ = index + Y;
                    if (++*(index + Z) == 0) *stackTopInner++ = index + Z;

                    counts[filterIdx] += CountExtensionsInner(
                        FilterDepth - 1, stackTopInner, stackTop2);

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
        }

        // Below FilterDepth: standard single-filter recursion with inline Work4
        Num CountExtensionsInner(int depth, byte** stackTop1, byte** stackTop2) {
            Num count = 0;
            byte** stackTopOriginal = stackTop1;
            while (stackTop1 != refStack) {
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
                    count += CountExtensionsInner(depth - 1, stackTopInner, stackTop2);
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
}

