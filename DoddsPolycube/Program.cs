using System.Text;

namespace DoddsPolycube;

// from https://oeis.org/A000162/a000162.cs.txt

// 59,795,121,480 free polycubes with 16 cells
// Runtime: 1:00.4 + 1:32.3 = 2:32.7
// parallel non-trivial and fix parallel trivial
// change order of adjacencyCounts
// de-if adjacencyCounts init
// Runtime: 28.2 + 1:32.2 = 2:00.4

// n                          a(n)
// 1                             1
// 2                             1
// 3                             2
// 4                             8
// 5                            29
// 6                           166
// 7                         1,023
// 8                         6,922
// 9                        48,311
// 10	                   346,543
// 11		             2,522,522
// 12		            18,598,427  ##  x24 = 446,362,248 (int limit is 2.1e9)
// 13		           138,462,649  ##  x24 = 3,322,710,384 (uint limit is 4.3e9)
// 14		         1,039,496,297 ## 1.0+0.9 = 1.9s
// 15		         7,859,514,470 ## 1.8+5.7 = 7.5s
// 16		        59,795,121,480 ## 6.1+38.8 = 44.9s
// 17		       457,409,613,979 ## 15+4:51 = 5:06 
// 18		     3,516,009,200,564 ## 55+37:36 = 38:31
// 19		    27,144,143,923,583 ## 2:23+04:59:14 = 05:01:38
// 20		   210,375,361,379,518 ## 8:50+1d 15:10:45=1d 15:19:35
// 21		 1,636,229,771,639,924 ## est. 12.7d
// 22		12,766,882,202,755,783 ## est. 100d  ##  x24 = 312e15 (ulong limit is 18.44e18)
// 23 est   99.9027e15 ## est. 2y  ##  x24 = 2.39766e18 (ulong is big enough)
// 24 est  783.7570e15 ## est. 16y ##  x24 = 18.8102e18 (UInt128 limit is 3.4e36, ulong is not quite big enough)
// 45 est    4.9048e36  ##  x24 = 117.715e36 (UInt128 is big enough)
// 46 est   37.5306e36  ##  x24 = 900.735e36 (need to switch to BigInteger)

using Num = ulong; // int, uint when n=13, ulong when n>=14, UInt128 when n>=24, BitInteger when n>45

using System.Diagnostics;
// using System.Numerics; // for BigInteger

public static class Program {
    private const int N = 16; // number of polycube cells. Need N >= 6 and > FilterDepth 

    // make sure type Num is big enough for a(N) * 24
    private const int FilterDepth = 5; // >=5 && <N

    // this is the maximum left stack length, which is the value being filtered to separate work
    private const int MaxLeftStackLen = 4 * (N - FilterDepth) - 2;

    private static int Main(string[] args) {
        if (N < 6) throw new InvalidOperationException("N must be at least 6");
        if (FilterDepth >= N) throw new InvalidOperationException("FilterDepth must be less than N");
        if (FilterDepth < 5) throw new InvalidOperationException("FilterDepth must be at least 5 for multithreading");
        string numType = typeof(Num).Name;
        bool numTypeOk = numType switch {
            "Int32" => N < 13,
            "UInt32" => N <= 13,
            "Int64" or "UInt64" => N <= 23,
            "Int128" or "UInt128" => N <= 45,
            "BigInteger" => true,
            _ => throw new InvalidOperationException($"Unknown Num type ({numType})")
        };
        if (!numTypeOk) throw new InvalidOperationException($"Num type ({numType}) is not big enough for N ({N})");

        /* n=16
           1,504,619 polycubes fixed under each orthogonal order 2 rotation - *3 = 4,513,857
           277 polycubes fixed under each orthogonal order 4 rotation - *6 = 1,662
           654,885 polycubes fixed under each short diagonal order 2 rotation - *6 = 3,929,310
           1,992 polycubes fixed under each long diagonal order 3 rotation - *8 = 15,936
           total count for nontrivial symmetries is 8,460,765 for polycubes with 16 cells - Elapsed: 00:01:56.8297006
         */

        HashSet<int> include = [..args.Select(int.Parse)];
        if (args.Length == 0)
            Console.WriteLine($"You can specify which symmetries to include as arguments (0-{
                MaxLeftStackLen}), (-1) for nontrivials");

        // enumerate the sum over the order 24 group of the size of the fix of each group element, and divide by 24
        // (Burnside's lemma)
        Num totalCount = 0;
        var swTotal = Stopwatch.StartNew();
        if (include.Count == 0 || include.Contains(-1)) {
            Console.WriteLine("Phase 1/2: started nontrivial symmetries");

            var filename = $"nontrivial_{N}.txt";
            if (File.Exists(filename)) {
                var lines = File.ReadAllLines(filename);
                var line0 = lines[0].Split(' ');
                totalCount = Num.Parse(line0[0]);
                //var elapsed = TimeSpan.Parse(line0[1]);
                for (int i = 1; i < lines.Length; i++)
                    Console.WriteLine(lines[i]);
            } else {
                Stopwatch sw = Stopwatch.StartNew();
                string[] descriptions = [
                    "orthogonal order 2",
                    "orthogonal order 4",
                    "short diagonal order 2",
                    "long diagonal order 3"
                ];
                int[] autClassSizes = [3, 6, 6, 8];
                int[][] matrixReps = [
                        [-1, 0, 0, 0, -1, 0, 0, 0, 1],
                        [0, -1, 0, 1, 0, 0, 0, 0, 1],
                        [0, 1, 0, 1, 0, 0, 0, 0, -1],
                        [0, 0, 1, 1, 0, 0, 0, 1, 0]
                    ],
                    affine1 = [[1, 0, 0], [1, 0, 0], [1, -1, 0], [1, 0, -1]],
                    affine2 = [[0, 1, 0], [0, 1, 0], [0, 0, 1], [0, 1, -1]],
                    biases = [[2 * N, 2 * N, 0], [2 * N, 0, 0], [0, 0, 2], [N - 1, 0, 1 - N]];
                List<Action> tasks = [];

                int completed = 0;
                var subCounts = new Num[4];
                var elapsed = TimeSpan.Zero;
                for (int sym = 0; sym < 4; sym++) {
                    int[] a1 = affine1[sym], a2 = affine2[sym], b = biases[sym];
                    for (int i = 1 - N; i <= N - 1; i++) {
                        for (int j = 1 - N + Math.Abs(i); j <= N - 1 - Math.Abs(i); j++) {
                            int[] matrixRep = matrixReps[sym],
                                affineShift = [
                                    i * a1[0] + j * a2[0] + b[0],
                                    i * a1[1] + j * a2[1] + b[1],
                                    i * a1[2] + j * a2[2] + b[2]
                                ];
                            int symCopy = sym, iCopy = i, jCopy = j; // copy, since lambda expression captures variable
                            tasks.Add(() => {
                                var sw = Stopwatch.StartNew();
                                var count = CountSymmetricPolycubes(matrixRep, affineShift);
                                sw.Stop();
                                lock (tasks) {
                                    completed++;
                                    subCounts[symCopy] += count;
                                    elapsed += sw.Elapsed;
                                }
                                Console.Write($"{completed}/{tasks.Count} {symCopy},{iCopy},{jCopy}     \r");
                            });
                        }
                    }
                }

                StringBuilder sb = new();
                var s = $"Starting {tasks.Count} tasks - start time: {DateTime.Now}";
                sb.AppendLine(s);
                Console.WriteLine(s);

                Parallel.Invoke(tasks.ToArray());
                // foreach (var task in tasks) task();

                for (int sym = 0; sym < 4; sym++) {
                    Num subCount = subCounts[sym], subCountMul = subCount * (Num)autClassSizes[sym];
                    s = $"    {sym}: {subCount:N0} polycubes fixed under each {descriptions[sym]} rotation - *{
                        autClassSizes[sym]} = {subCountMul:N0}";
                    sb.AppendLine(s);
                    Console.WriteLine(s);

                    totalCount += subCountMul;
                }

                s = $"total count for nontrivial symmetries is {totalCount:N0} for polycubes with {
                    N} cells - Elapsed: {sw.Elapsed}, CPU time: {elapsed}";
                sb.AppendLine(s);
                Console.WriteLine(s);

                sb.Insert(0, $"{totalCount} {sw.Elapsed}{Environment.NewLine}");
                File.WriteAllText(filename, sb.ToString());
            }
        }
        Console.WriteLine();
        {
            Console.WriteLine("Phase 2/2: started trivial symmetries");
            var sw = Stopwatch.StartNew();
            List<Action> tasks = [];
            Num subCount = 0;
            var cpuTime = TimeSpan.Zero;
            for (int j = 0, completed = 0; j <= MaxLeftStackLen; j++) {
                if (include.Count != 0 && !include.Contains(j)) continue;
                var filename = $"trivial_{N}_{MaxLeftStackLen}_{j}.txt";
                if (File.Exists(filename)) {
                    var lines = File.ReadAllLines(filename);
                    var line0 = lines[0].Split(' ');
                    subCount += Num.Parse(line0[0]);
                    cpuTime += TimeSpan.Parse(line0[1]);
                    for (int i = 1; i < lines.Length; i++)
                        Console.WriteLine(lines[i]);
                } else {
                    int filter = j; // copy, since lambda expression captures the variable
                    tasks.Add(() => {
                        var sw = Stopwatch.StartNew();
                        var count = CountExtensionsSubsetUnsafe(filter);
                        sw.Stop();
                        lock (tasks) {
                            completed++;
                            subCount += count;
                            cpuTime += sw.Elapsed;
                        }
                        var s = $"[{completed}/{tasks.Count}] #{filter} count={count:N0} elapsed={sw.Elapsed}";
                        File.WriteAllText(filename, $"{count} {sw.Elapsed}{Environment.NewLine}" + s);
                        Console.WriteLine(s);
                    });
                }
            }

            Console.WriteLine($"Starting {tasks.Count} tasks - start time: {DateTime.Now}");

            Parallel.Invoke(tasks.ToArray());
            // foreach (var task in tasks) task();

            Console.WriteLine($"{subCount:N0} polycubes with {
                N} cells (number of polycubes fixed by trivial symmetry) - Elapsed: {sw.Elapsed}, CPU time: {cpuTime}");

            totalCount += subCount;
        }

        if (include.Count == 0) {
            Console.WriteLine();
            totalCount /= 24;
            Console.WriteLine($"{totalCount:N0} free polycubes with {N} cells - Elapsed: {swTotal.Elapsed}");
        }
        return 0;
    }

    private static Num CountSymmetricPolycubes(int[] linearMap, int[] affineShift) {
        const int XMul = 1, YMul = 2 * N + 1, ZMul = YMul * (2 * N + 1), Size = (N + 2) * ZMul;
        var adjacencyCounts = new byte[Size];
        Array.Fill(adjacencyCounts, (byte)1, 0, ZMul + N * YMul + N + 1);

        HashSet<(int, int, int)> requiredCells = [];
        Stack<(int, int, int)> recoveryStack = new(), extensionStack = new();
        extensionStack.Push((N, N, 1));
        return CountExtensions(N);

        Num CountExtensions(int cellsToAdd) {
            cellsToAdd--;
            Num count = 0;
            int originalLength = extensionStack.Count;
            while (extensionStack.Count > 0) {
                int x, y, z;
                recoveryStack.Push((x, y, z) = extensionStack.Pop());
                int xyz = x * XMul + y * YMul + z * ZMul;

                bool existingRequirement = requiredCells.Remove((x, y, z));
                if (!existingRequirement) {
                    if (cellsToAdd < requiredCells.Count)
                        continue; // number of required cells will only grow, so if already impossible, go to next
                    for (int tempX = x, tempY = y, tempZ = z;;) { // works for general transformations of finite order
                        (tempX, tempY, tempZ) = (
                            linearMap[0] * tempX + linearMap[1] * tempY + linearMap[2] * tempZ + affineShift[0],
                            linearMap[3] * tempX + linearMap[4] * tempY + linearMap[5] * tempZ + affineShift[1],
                            linearMap[6] * tempX + linearMap[7] * tempY + linearMap[8] * tempZ + affineShift[2]);
                        if (x == tempX && y == tempY && z == tempZ) break;
                        requiredCells.Add((tempX, tempY, tempZ));
                    }
                }

                // if there are too many required cells, then no valid polycubes are possible
                if (cellsToAdd >= requiredCells.Count) {
                    if (cellsToAdd == 0) count++;
                    else {
                        int innerOriginalLength = extensionStack.Count;
                        if (adjacencyCounts[xyz - XMul]++ == 0) extensionStack.Push((x - 1, y, z));
                        if (adjacencyCounts[xyz - YMul]++ == 0) extensionStack.Push((x, y - 1, z));
                        if (adjacencyCounts[xyz - ZMul]++ == 0) extensionStack.Push((x, y, z - 1));
                        if (adjacencyCounts[xyz + XMul]++ == 0) extensionStack.Push((x + 1, y, z));
                        if (adjacencyCounts[xyz + YMul]++ == 0) extensionStack.Push((x, y + 1, z));
                        if (adjacencyCounts[xyz + ZMul]++ == 0) extensionStack.Push((x, y, z + 1));

                        count += CountExtensions(cellsToAdd);

                        --adjacencyCounts[xyz - XMul];
                        --adjacencyCounts[xyz - YMul];
                        --adjacencyCounts[xyz - ZMul];
                        --adjacencyCounts[xyz + XMul];
                        --adjacencyCounts[xyz + YMul];
                        --adjacencyCounts[xyz + ZMul];
                        while (extensionStack.Count != innerOriginalLength)
                            extensionStack.Pop(); // maybe replace this w/ custom stack to avoid this loop
                    }
                }

                if (existingRequirement) {
                    requiredCells.Add((x, y, z));
                    break; // this required cell will no longer be available in the extension stack,
                    // so no more valid polycubes are possible in this branch
                }

                for (int tempX = x, tempY = y, tempZ = z;;) {
                    (tempX, tempY, tempZ) = (
                        linearMap[0] * tempX + linearMap[1] * tempY + linearMap[2] * tempZ + affineShift[0],
                        linearMap[3] * tempX + linearMap[4] * tempY + linearMap[5] * tempZ + affineShift[1],
                        linearMap[6] * tempX + linearMap[7] * tempY + linearMap[8] * tempZ + affineShift[2]);
                    if (x == tempX && y == tempY && z == tempZ) break;
                    requiredCells.Remove((tempX, tempY, tempZ));
                }
            }
            while (extensionStack.Count != originalLength)
                extensionStack.Push(recoveryStack.Pop());
            return count;
        }
    }

    private static Num CountExtensionsSubsetUnsafe(int filter) {
        unsafe {
            // set X<Y<Z such that aX+bY+cZ = 0 implies a = b = c = 0 or |a|+|b|+|c| > n
            const int X = (N + 5) / 4 * ((N + 5) / 4 * 3 - 2);
            // trivial choice is X = 1, Y = n, Z = n * n.
            // A simple reduction is X = 1, Y = n, Z = n * (n / 2) + (n + 1) / 2
            const int Y = X + 1;
            // minimising Z is memory efficient. Unclear if this noticably affects performance.
            // Z ~ 3/16 n^2 is the best I can find for arbitrary n
            const int Z = X + (N + 5) / 4 * 3;
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
                    byte** stackTopInner = stackTop1;

                    if (++*(index - X) == 0) *stackTopInner++ = index - X;
                    if (++*(index - Y) == 0) *stackTopInner++ = index - Y;
                    if (++*(index - Z) == 0) *stackTopInner++ = index - Z;
                    if (++*(index + X) == 0) *stackTopInner++ = index + X;
                    if (++*(index + Y) == 0) *stackTopInner++ = index + Y;
                    if (++*(index + Z) == 0) *stackTopInner++ = index + Z;

                    if (depth == 4) {
                        byte** stackTop = stackTopInner;
                        const int X2 = X + X,
                            Y2 = Y + Y,
                            Z2 = Z + Z,
                            sYX = Y + X,
                            sZX = Z + X,
                            sZY = Z + Y,
                            dYX = Y - X,
                            dZX = Z - X,
                            dZY = Z - Y;
                        int length = (int)(stackTop - refStack);
                        count += (Num)(length * (length - 1) * (length - 2) / 6);
                        byte** stackTopTemp = stackTop;
                        for (int lengthPlus = (length << 1) - 511; stackTopTemp != refStack;) {
                            byte* i = *--stackTopTemp;
                            int neighbours = 0, subCount = 128;
                            if (*(i - X) > 127) {
                                count += --*(i - X);
                                neighbours++;
                                subCount += *(i - X2) + *(i - sYX) + *(i - sZX) + *(i + dYX) + *(i + dZX);
                            }
                            if (*(i - Y) > 127) {
                                count += --*(i - Y);
                                neighbours++;
                                subCount += *(i - Y2) + *(i - sYX) + *(i - sZY) + *(i - dYX) + *(i + dZY);
                            }
                            if (*(i - Z) > 127) {
                                count += --*(i - Z);
                                neighbours++;
                                subCount += *(i - Z2) + *(i - sZX) + *(i - sZY) + *(i - dZX) + *(i - dZY);
                            }
                            if (*(i + X) > 127) {
                                count += --*(i + X);
                                neighbours++;
                                subCount += *(i + X2) + *(i + sYX) + *(i + sZX) + *(i - dYX) + *(i - dZX);
                            }
                            if (*(i + Y) > 127) {
                                count += --*(i + Y);
                                neighbours++;
                                subCount += *(i + Y2) + *(i + sYX) + *(i + sZY) + *(i + dYX) + *(i - dZY);
                            }
                            if (*(i + Z) > 127) {
                                count += --*(i + Z);
                                neighbours++;
                                subCount += *(i + Z2) + *(i + sZX) + *(i + sZY) + *(i + dZX) + *(i + dZY);
                            }
                            count += (Num)((subCount >> 8) + (neighbours * (neighbours + lengthPlus) >> 1));
                        }
                        while (stackTop != refStack) {
                            byte* i = *--stackTop;
                            *(i - X) |= (byte)(*(i - X) >> 4);
                            *(i - Y) |= (byte)(*(i - Y) >> 4);
                            *(i - Z) |= (byte)(*(i - Z) >> 4);
                            *(i + X) |= (byte)(*(i + X) >> 4);
                            *(i + Y) |= (byte)(*(i + Y) >> 4);
                            *(i + Z) |= (byte)(*(i + Z) >> 4);
                        }
                    } else if (depth != FilterDepth || stackTop1 - refStack == filter)
                        // if multithreading is not wanted, remove "if (condition)" from this else statement
                        count += CountExtensions(depth - 1, stackTopInner, stackTop2);

                    --*(index - X);
                    --*(index - Y);
                    --*(index - Z);
                    --*(index + X);
                    --*(index + Y);
                    --*(index + Z);
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

    private static Num CountExtensionsSubset(int filter) {
        // set X<Y<Z such that aX+bY+cZ = 0 implies a = b = c = 0 or |a|+|b|+|c| > n
        const int X = (N + 5) / 4 * ((N + 5) / 4 * 3 - 2);
        // trivial choice is X = 1, Y = n, Z = n * n.
        // A simple reduction is X = 1, Y = n, Z = n * (n / 2) + (n + 1) / 2
        const int Y = X + 1;
        // minimising Z is memory efficient. Unclear if this noticably affects performance.
        // Z ~ 3/16 n^2 is the best I can find for arbitrary n
        const int Z = X + (N + 5) / 4 * 3;
        // could use ints or shorts as offsets to save memory, but it's faster to directly store the
        // pointers to avoid adding pointer offsets at every lookup
        byte[] byteBoard = new byte[(N + 2) * Z];
        // the first Z + 1 bytes are disallowed extensions; first Z are less than the minimum,
        // last 1 due to edge case of initial polycube having no neighbours
        Array.Fill(byteBoard, (byte)255, Z + 1, byteBoard.Length - Z - 1);

        // total length of the two stacks is at most 4n-9. One stack grows from the left, the other
        // stack grows from the right
        int[] refStack = new int[(N - 2) * 4];
        // seeded with first index of the byte board as the only allowed extension
        refStack[0] = Z;

        return CountExtensions(N, 1, refStack.Length);

        Num CountExtensions(int depth, int stackPtr, int stackLimit) {
            Num count = 0;
            int stackTopOriginal = stackPtr;
            while (stackPtr != 0) {
                int index = refStack[--stackPtr];
                int stackTopInner = stackPtr;

                if (++byteBoard[index - X] == 0) refStack[stackTopInner++] = index - X;
                if (++byteBoard[index - Y] == 0) refStack[stackTopInner++] = index - Y;
                if (++byteBoard[index - Z] == 0) refStack[stackTopInner++] = index - Z;
                if (++byteBoard[index + X] == 0) refStack[stackTopInner++] = index + X;
                if (++byteBoard[index + Y] == 0) refStack[stackTopInner++] = index + Y;
                if (++byteBoard[index + Z] == 0) refStack[stackTopInner++] = index + Z;

                if (depth == 4) {
                    int stackTop = stackTopInner;
                    const int X2 = X + X, Y2 = Y + Y, Z2 = Z + Z;
                    const int sYX = Y + X, sZX = Z + X, sZY = Z + Y;
                    const int dYX = Y - X, dZX = Z - X, dZY = Z - Y;
                    int length = stackTop;
                    count += (Num)(length * (length - 1) * (length - 2) / 6);
                    int stackTopTemp = stackTop;
                    for (int lengthPlus = (length << 1) - 511; stackTopTemp != 0;) {
                        int i = refStack[--stackTopTemp];
                        int neighbours = 0, subCount = 128;
                        if (byteBoard[i - X] > 127) {
                            count += --byteBoard[i - X];
                            neighbours++;
                            subCount += byteBoard[i - X2] +
                                        byteBoard[i - sYX] +
                                        byteBoard[i - sZX] +
                                        byteBoard[i + dYX] +
                                        byteBoard[i + dZX];
                        }
                        if (byteBoard[i - Y] > 127) {
                            count += --byteBoard[i - Y];
                            neighbours++;
                            subCount += byteBoard[i - Y2] +
                                        byteBoard[i - sYX] +
                                        byteBoard[i - sZY] +
                                        byteBoard[i - dYX] +
                                        byteBoard[i + dZY];
                        }
                        if (byteBoard[i - Z] > 127) {
                            count += --byteBoard[i - Z];
                            neighbours++;
                            subCount += byteBoard[i - Z2] +
                                        byteBoard[i - sZX] +
                                        byteBoard[i - sZY] +
                                        byteBoard[i - dZX] +
                                        byteBoard[i - dZY];
                        }
                        if (byteBoard[i + X] > 127) {
                            count += --byteBoard[i + X];
                            neighbours++;
                            subCount += byteBoard[i + X2] +
                                        byteBoard[i + sYX] +
                                        byteBoard[i + sZX] +
                                        byteBoard[i - dYX] +
                                        byteBoard[i - dZX];
                        }
                        if (byteBoard[i + Y] > 127) {
                            count += --byteBoard[i + Y];
                            neighbours++;
                            subCount += byteBoard[i + Y2] +
                                        byteBoard[i + sYX] +
                                        byteBoard[i + sZY] +
                                        byteBoard[i + dYX] +
                                        byteBoard[i - dZY];
                        }
                        if (byteBoard[i + Z] > 127) {
                            count += --byteBoard[i + Z];
                            neighbours++;
                            subCount += byteBoard[i + Z2] +
                                        byteBoard[i + sZX] +
                                        byteBoard[i + sZY] +
                                        byteBoard[i + dZX] +
                                        byteBoard[i + dZY];
                        }
                        count += (Num)((subCount >> 8) + (neighbours * (neighbours + lengthPlus) >> 1));
                    }
                    while (stackTop != 0) {
                        int i = refStack[--stackTop];
                        byteBoard[i - X] |= (byte)(byteBoard[i - X] >> 4);
                        byteBoard[i - Y] |= (byte)(byteBoard[i - Y] >> 4);
                        byteBoard[i - Z] |= (byte)(byteBoard[i - Z] >> 4);
                        byteBoard[i + X] |= (byte)(byteBoard[i + X] >> 4);
                        byteBoard[i + Y] |= (byte)(byteBoard[i + Y] >> 4);
                        byteBoard[i + Z] |= (byte)(byteBoard[i + Z] >> 4);
                    }
                } else if (depth != FilterDepth || stackPtr == filter) {
                    // if multithreading is not wanted, remove "if (condition)" from this else statement
                    count += CountExtensions(depth - 1, stackTopInner, stackLimit);
                }

                --byteBoard[index - X];
                --byteBoard[index - Y];
                --byteBoard[index - Z];
                --byteBoard[index + X];
                --byteBoard[index + Y];
                --byteBoard[index + Z];

                // doing this push before the recursion would add one extra unnecessary element to the stack
                // at each level of recursion
                refStack[--stackLimit] = index;
            }

            while (stackPtr != stackTopOriginal) {
                refStack[stackPtr++] = refStack[stackLimit++];
            }
            return count;
        }
    }
}

// written by Stanley Dodds who says:
//   I can probably supply a full proof of correctness if needed, but I have nothing written down right now
//   so currently, source: dude, just trust me (or try reading my code)
//   Of course, you can compare up to n = 19 with previous results and see that it's definitely doing something right