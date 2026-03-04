namespace DoddsPolycube;

using Num = ulong;

public partial class Program {
    /// <summary>
    /// Prints a table analyzing the tradeoffs of different FilterDepth values for the current N.
    /// Helps determine the optimal FilterDepth before a long computation.
    /// </summary>
    private static void PrintFilterDepthAnalysis(int threads) {
        Console.WriteLine($"FilterDepth analysis for N={N}, {threads} threads");
        Console.WriteLine($"Current setting: FilterDepth={FilterDepth}, MaxLeftStackLen={MaxLeftStackLen}");
        Console.WriteLine();
        Console.WriteLine("FD  Filters  Per-Thread  Inner  Outer  Notes");
        Console.WriteLine("--  -------  ----------  -----  -----  -----");

        for (int fd = 5; fd < N; fd++) {
            int maxLeftStack = 4 * (N - fd) - 2;
            int filters = maxLeftStack + 1;
            int innerDepth = fd - 5;
            int outerDepth = N - fd - 1;
            double filtersPerThread = (double)filters / threads;

            string notes = "";
            if (filters < threads)
                notes += $"⚠ idle threads ({filters}<{threads})  ";
            if (filtersPerThread < 2)
                notes += "⚠ poor load balance  ";
            if (outerDepth <= 0)
                notes += "⚠ no checkpointing  ";
            if (fd == FilterDepth)
                notes += "◄ current  ";

            Console.WriteLine($"{fd,2}  {filters,7}  {filtersPerThread,10:F1}  {innerDepth,5}  {outerDepth,5}  {notes}");
        }

        Console.WriteLine();
        Console.WriteLine("To find the optimal FD empirically:");
        Console.WriteLine("  1. Set N to a small value (16-18) that completes in minutes");
        Console.WriteLine("  2. Run: --fdbenchmark -P  (tests all FD values automatically)");
        Console.WriteLine("  3. Use the winner for your production N=23 run");
        Console.WriteLine();
        Console.WriteLine("From N=16 benchmarks (8 threads, M2 Max), the sweet spot is FD=6-10.");
        Console.WriteLine("The plateau is broad — any value in range gives ~13% speedup over FD=5.");
        Console.WriteLine("Lower FD = more filters = better load balance + finer checkpoints.");
        Console.WriteLine("Higher FD = more inner depth = slightly faster per-operation.");
        Console.WriteLine("Recommendation for N=23 with 8 threads: FD=7 or FD=8.");
    }

    /// <summary>
    /// Benchmarks all valid FD values by running the full computation at the current N.
    /// For small N (16-18), this completes in minutes. For large N, use --fdquick instead.
    /// </summary>
    private static void RunFilterDepthBenchmark(int threads, bool quickMode) {
        Console.WriteLine($"FD benchmark: N={N}, {threads} threads{(quickMode ? " (quick mode: filter 0 only, single-threaded)" : "")}");
        int maxFd = Math.Min(N - 2, 15);
        Console.WriteLine($"Testing FilterDepth values from 5 to {maxFd}...");
        Console.WriteLine();

        int boardSize = (N + 2) * Z;
        int stackSize = (N - 2) * 4;

        string bestResult = "";
        double bestTime = double.MaxValue;

        for (int fd = 5; fd <= maxFd; fd++) {
            int maxLeft = 4 * (N - fd) - 2;
            int filterCount = maxLeft + 1;
            if (!quickMode && filterCount < threads) {
                Console.WriteLine($"  FD={fd,2}: {filterCount} filters < {threads} threads — skipping");
                continue;
            }

            if (quickMode) {
                // Quick mode: run filter 0 only, single-threaded
                var board = new byte[boardSize];
                var stack = new int[stackSize];
                Array.Fill(board, (byte)255, Z + 1, boardSize - Z - 1);
                stack[0] = Z;
                var counts = new Num[filterCount];
                var filterSet = new HashSet<int> { 0 };

                var sw = System.Diagnostics.Stopwatch.StartNew();
                RunOuterLoop(board, stack, stackSize, filterSet, counts, fd);
                sw.Stop();

                double secs = sw.Elapsed.TotalSeconds;
                string marker = secs < bestTime ? " ◄ best so far" : "";
                if (secs < bestTime) { bestTime = secs; bestResult = $"FD={fd}"; }
                Console.WriteLine($"  FD={fd,2}: {sw.Elapsed:mm\\:ss\\.f}  (filter 0 only)  count={counts[0]:N0}{marker}");
            } else {
                // Full mode: run all filters with parallelism
                int actualPaths = Math.Min(threads, filterCount);
                var groups = new HashSet<int>[actualPaths];
                for (int p = 0; p < actualPaths; p++) groups[p] = [];
                for (int f = 0; f < filterCount; f++)
                    groups[f % actualPaths].Add(f);

                var allCounts = new Num[filterCount];

                var sw = System.Diagnostics.Stopwatch.StartNew();
                Parallel.Invoke(groups.Select(fs => (Action)(() => {
                    var board = new byte[boardSize];
                    var stack = new int[stackSize];
                    Array.Fill(board, (byte)255, Z + 1, boardSize - Z - 1);
                    stack[0] = Z;
                    var counts = new Num[filterCount];
                    RunOuterLoop(board, stack, stackSize, fs, counts, fd);
                    lock (allCounts) {
                        for (int f = 0; f < filterCount; f++)
                            allCounts[f] += counts[f];
                    }
                })).ToArray());
                sw.Stop();

                Num total = 0;
                for (int f = 0; f < filterCount; f++) total += allCounts[f];

                double secs = sw.Elapsed.TotalSeconds;
                string marker = secs < bestTime ? " ◄ best so far" : "";
                if (secs < bestTime) { bestTime = secs; bestResult = $"FD={fd}"; }
                Console.WriteLine($"  FD={fd,2}: {sw.Elapsed:mm\\:ss\\.f}  ({filterCount} filters, {filterCount * 1.0 / threads:F1}/thread)  count={total:N0}{marker}");
            }
        }

        Console.WriteLine();
        Console.WriteLine($"Optimal: {bestResult} at {bestTime:F1}s");
        Console.WriteLine($"Set FilterDepth to this value in Program.cs and recompile for production runs.");
    }

    /// <summary>
    /// Runs the managed outer loop + unsafe inner dispatch.
    /// Takes filterDepth as a runtime parameter for benchmarking different values.
    /// </summary>
    private static unsafe void RunOuterLoop(byte[] byteBoard, int[] refStack,
        int stackSize, HashSet<int> filterSet, Num[] filterCounts, int filterDepth) {

        int[] callStack = new int[stackSize];
        int depth = N, stackPtr = 1, stackTopOriginal = 1, stackLimit = stackSize;
        int index = 0, callStackPtr = 0;
        bool popping = false;

        // Pre-allocate pointer stack once, reused for each dispatch
        fixed (byte* boardBase = byteBoard) {
        byte** ptrStack = stackalloc byte*[stackSize];

        for (;;) {
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

                if (depth == filterDepth + 1) {
                    // Dispatch using pre-allocated pointer stack
                    for (int i = 0; i < stackTopInner; i++)
                        ptrStack[i] = boardBase + refStack[i];

                    byte** ptrBase = ptrStack;
                    byte** sTop1 = ptrStack + stackTopInner;
                    byte** sTop2 = ptrStack + stackLimit;

                    byte** sTopOrig = sTop1;
                    while (sTop1 != ptrBase) {
                        byte* idx = *--sTop1;
                        int filterIdx = (int)(sTop1 - ptrBase);
                        if (filterSet.Contains(filterIdx)) {
                            byte** sTopInner = sTop1;
                            if (++*(idx - Z) == 0) *sTopInner++ = idx - Z;
                            if (++*(idx - Y) == 0) *sTopInner++ = idx - Y;
                            if (++*(idx - X) == 0) *sTopInner++ = idx - X;
                            if (++*(idx + X) == 0) *sTopInner++ = idx + X;
                            if (++*(idx + Y) == 0) *sTopInner++ = idx + Y;
                            if (++*(idx + Z) == 0) *sTopInner++ = idx + Z;

                            filterCounts[filterIdx] += CountExtInner(
                                ptrBase, filterDepth - 1, sTopInner, sTop2);

                            --*(idx - Z); --*(idx - Y); --*(idx - X);
                            --*(idx + X); --*(idx + Y); --*(idx + Z);
                        }
                        *--sTop2 = idx;
                    }
                    while (sTop1 != sTopOrig)
                        *sTop1++ = *sTop2++;
                    for (int i = 0; i < stackTopInner; i++)
                        refStack[i] = (int)(ptrStack[i] - boardBase);
                } else if (depth > filterDepth + 1) {
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
        } // end fixed
    }
}





