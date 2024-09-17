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
// 14		         1,039,496,297
// 15		         7,859,514,470
// 16		        59,795,121,480 ## 6.1+38.8 = 44.9
// 17		       457,409,613,979 ## 15.6+4:51.4 = 5:07.0 
// 18		     3,516,009,200,564 ## 55.2+37:35.9 = 38:31.0
// 19		    27,144,143,923,583 ## 2:23.4+04:59:14.4 = 05:01:37.8
// 20		   210,375,361,379,518
// 21		 1,636,229,771,639,924
// 22		12,766,882,202,755,783  ##  x24 = 312e15 (ulong limit is 18.44e18)
// 23 est   99.9027e15  ##  x24 = 2.39766e18 (ulong is big enough)
// 24 est  783.7570e15  ##  x24 = 18.8102e18 (UInt128 limit is 3.4e36, ulong is not quite big enough)
// 45 est    4.9048e36  ##  x24 = 117.715e36 (UInt128 is big enough)
// 46 est   37.5306e36  ##  x24 = 900.735e36 (need to switch to BigInteger)

using Num = long; // int, uint when n=13, ulong when n>=14, UInt128 when n>=24, BitInteger when n>45

using System.Diagnostics;
// using System.Numerics; // for BigInteger

public static class Program {
    private const int N = 16; // number of polycube cells. Need N >= 6 and > FilterDepth 
    // make sure type Num is big enough for a(N) * 24
    private const int FilterDepth = 5; // >=5 && <N

    private static void Main() {
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
            _ => false
        };
        if (!numTypeOk) throw new InvalidOperationException($"Num type ({numType}) is not big enough for N ({N})");

        /* n=16
           1,504,619 polycubes fixed under each orthogonal order 2 rotation - *3 = 4,513,857
           277 polycubes fixed under each orthogonal order 4 rotation - *6 = 1,662
           654,885 polycubes fixed under each short diagonal order 2 rotation - *6 = 3,929,310
           1,992 polycubes fixed under each long diagonal order 3 rotation - *8 = 15,936
           total count for nontrivial symmetries is 8,460,765 for polycubes with 16 cells - Elapsed: 00:01:56.8297006
         */

        // enumerate the sum over the order 24 group of the size of the fix of each group element, and divide by 24
        // (Burnside's lemma)
        Num totalCount = 0;
        Stopwatch swTotal = Stopwatch.StartNew();
        {
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
            ], affine1 = [[1, 0, 0], [1, 0, 0], [1, -1, 0], [1, 0, -1]],
                affine2 = [[0, 1, 0], [0, 1, 0], [0, 0, 1], [0, 1, -1]],
                biases = [[2 * N, 2 * N, 0], [2 * N, 0, 0], [0, 0, 2], [N - 1, 0, 1 - N]];
            List<Task<Num>> tasks = [];
            int spindex = 0;
            const string spinner = @"|/-\";
            
            for (int sym = 0; sym < 4; sym++) {
                int[] a1 = affine1[sym], a2 = affine2[sym], b = biases[sym];
                for (int i = 1 - N; i <= N - 1; i++) {
                    for (int j = 1 - N + Math.Abs(i); j <= N - 1 - Math.Abs(i); j++) {
                        var matrixRep = matrixReps[sym];
                        var affineShift = new[] {
                            i * a1[0] + j * a2[0] + b[0],
                            i * a1[1] + j * a2[1] + b[1],
                            i * a1[2] + j * a2[2] + b[2]};
                        tasks.Add(Task<Num>.Factory.StartNew(() => {
                            Num count = CountSymmetricPolycubes(matrixRep, affineShift);
                            Interlocked.Increment(ref spindex);
                            Console.Write($"{spinner[spindex % 4]}\b");
                            return count;
                        }));
                    }
                }
            }

            Console.WriteLine("Phase 1/2: started nontrivial symmetries - {0} tasks - start time: {1}",
                tasks.Count, DateTime.Now);

            Task.WhenAll(tasks).Wait();

            for (int sym = 0, taskIndex = 0; sym < 4; sym++) {
                Num subCount = 0;
                for (int i = 1 - N; i <= N - 1; i++)
                    for (int j = 1 - N + Math.Abs(i); j <= N - 1 - Math.Abs(i); j++)
                        subCount += tasks[taskIndex++].Result;
                
                Console.Write($"    {sym}: {subCount:N0} polycubes fixed under each {descriptions[sym]} rotation - ");
                subCount *= autClassSizes[sym];
                Console.WriteLine($"*{autClassSizes[sym]} = {subCount:N0}");
                totalCount += subCount;
            }
            
            Console.WriteLine("total: {0:N0} polycubes with {1} cells with nontrivial symmetries - Elapsed: {2}",
                totalCount, N, sw.Elapsed);
        }
        Console.WriteLine();
        {
            Stopwatch sw = Stopwatch.StartNew();
            // this is the maximum left stack length, which is the value being filtered to separate work
            int maxLeftStackLen = 4 * (N - FilterDepth) - 2;
            var tasks = new Task<Num>[maxLeftStackLen + 1];
            
            Console.WriteLine("Phase 2/2: started trivial symmetries - {0} tasks - start time: {1}",
                tasks.Length, DateTime.Now);
            
            Num subCount = 0;
            for (int j = 0, completed = 0; j < tasks.Length; j++) {
                int filter = maxLeftStackLen--; // copy, since lambda expression captures the variable
                tasks[j] = Task<Num>.Factory.StartNew(() => {
                    Num count = CountExtensionsSubset(filter);
                    lock (tasks) {
                        completed++;
                        subCount += count;
                    }
                    Console.Write($"[{completed}/{tasks.Length}] subcount={subCount:N0} elapsed={
                        sw.Elapsed}     \r");
                    return count;
                });
            }

            Task.WhenAll(tasks).Wait();

            Num subCount2 = 0;
            foreach (var t in tasks)
                subCount2 += t.Result;
            
            Console.WriteLine(
                "{0:N0} polycubes with {1} cells (number of polycubes fixed by trivial symmetry) - Elapsed: {2}",
                subCount2, N, sw.Elapsed);
            
            totalCount += subCount2;
        }
        
        Console.WriteLine();
        totalCount /= 24;
        Console.WriteLine("{0:N0} free polycubes with {1} cells - Elapsed: {2}", totalCount, N, swTotal.Elapsed);
    }

    private static Num CountSymmetricPolycubes(int[] linearMap, int[] affineShift) {
        var adjacencyCounts = new byte[2 * N + 1, 2 * N + 1, N + 2];
        int z = 0, y = 0, x;
        for (; y < 2 * N + 1; y++)
            for (x = 0; x < 2 * N + 1; x++)
                adjacencyCounts[x, y, z] = 1;
        for (z++, y = 0; y < N; y++)
            for (x = 0; x < 2 * N + 1; x++)
                adjacencyCounts[x, y, z] = 1;
        // y = n;
        for (x = 0; x <= N; x++)
            adjacencyCounts[x, y, z] = 1;
        
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
                        if (adjacencyCounts[x - 1, y, z]++ == 0) extensionStack.Push((x - 1, y, z));
                        if (adjacencyCounts[x, y - 1, z]++ == 0) extensionStack.Push((x, y - 1, z));
                        if (adjacencyCounts[x, y, z - 1]++ == 0) extensionStack.Push((x, y, z - 1));
                        if (adjacencyCounts[x + 1, y, z]++ == 0) extensionStack.Push((x + 1, y, z));
                        if (adjacencyCounts[x, y + 1, z]++ == 0) extensionStack.Push((x, y + 1, z));
                        if (adjacencyCounts[x, y, z + 1]++ == 0) extensionStack.Push((x, y, z + 1));

                        count += CountExtensions(cellsToAdd);

                        --adjacencyCounts[x - 1, y, z];
                        --adjacencyCounts[x, y - 1, z];
                        --adjacencyCounts[x, y, z - 1];
                        --adjacencyCounts[x + 1, y, z];
                        --adjacencyCounts[x, y + 1, z];
                        --adjacencyCounts[x, y, z + 1];
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

    private static Num CountExtensionsSubset(int filter) {
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
                        count += length * (length - 1) * (length - 2) / 6;
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
                            count += (subCount >> 8) + (neighbours * (neighbours + lengthPlus) >> 1);
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
                    }  else if (depth != FilterDepth || stackTop1 - refStack == filter)
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
}

// written by Stanley Dodds who says:
//   I can probably supply a full proof of correctness if needed, but I have nothing written down right now
//   so currently, source: dude, just trust me (or try reading my code)
//   Of course, you can compare up to n = 19 with previous results and see that it's definitely doing something right