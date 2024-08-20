namespace DoddsPolycube;

// from https://oeis.org/A000162/a000162.cs.txt

// 59,795,121,480 free polycubes with 16 cells
// Runtime: 1:00.4 + 1:32.3 = 2:32.7
// parallel non-trivial and fix parallel trivial
// change order of adjacencyCounts
// de-if adjencyCounts init
// Runtime: 28.2 + 1:32.2 = 2:00.4

using System.Diagnostics;

internal class Program
{
    private static int n = 23; //number of polycube cells. Need n >= 4 if single threading, or n >= filterDepth >= 5 if multithreading (I think)

    private static int Main(string[] args)
    {
        if (args.Length > 0)
        {
            n = int.Parse(args[0]);
        }
        if (n < 9)
        {
            Console.WriteLine("n must be at least 9");
            return 1;
        }

        const int filterDepth = 12;
        Console.WriteLine("Computing number of free polycubes with {0} cells using filter depth of {1}", n, filterDepth);

        long count = 0;//enumerate the sum over the order 24 group of the size of the fix of each group element, and divide by 24 (Burnside's lemma)
        TimeSpan totalTime = TimeSpan.Zero;
        {
            Console.WriteLine("Phase 1/2: Computing non-trivial symmetries");
            Stopwatch sw = Stopwatch.StartNew();
            string[] descriptions = ["orthogonal order 2 rotation", "orthogonal order 4 rotation", "short diagonal order 2 rotation", "long diagonal order 3 rotation"];
            int[] autClassSizes = [3, 6, 6, 8];
            int[][] matrixReps = [[-1, 0, 0, 0, -1, 0, 0, 0, 1], [0, -1, 0, 1, 0, 0, 0, 0, 1], [0, 1, 0, 1, 0, 0, 0, 0, -1], [0, 0, 1, 1, 0, 0, 0, 1, 0]];
            int[][] affine1 = [[1, 0, 0], [1, 0, 0], [1, -1, 0], [1, 0, -1]];
            int[][] affine2 = [[0, 1, 0], [0, 1, 0], [0, 0, 1], [0, 1, -1]];
            int[][] biases = [[2 * n, 2 * n, 0], [2 * n, 0, 0], [0, 0, 2], [n - 1, 0, 1 - n]];
            bool[,,] done = new bool[4, 2 * n, 2 * n];
            long[] subcounts = new long[4];
            int progress = 0, total = 4 * 2 * n * 2 * n;
            TimeSpan highestTime = TimeSpan.Zero;
            long subcount = 0;
            Enumerable.Range(0, 4).SelectMany(sym => Enumerable.Range(0, 2 * n).SelectMany(y => Enumerable.Range(0, 2 * n).Select(x => (sym, y, x)))).AsParallel().ForAll(e =>
            {
                int sym = e.sym, y = e.y, x = e.x;
                int i = 1 - n + y, j = 1 - n + Math.Abs(i) + x;
                string fileName = $"{n}-{sym}_{i}_{j}.txt";
                if (File.Exists(fileName))
                {
                    string[] fileParts = File.ReadAllText(fileName).Split(' ');
                    long sc = long.Parse(fileParts[0]);
                    if (fileParts.Length > 1 && TimeSpan.Parse(fileParts[1]) > highestTime) highestTime = TimeSpan.Parse(fileParts[1]);
                    Interlocked.Add(ref subcounts[sym], sc);
                    Interlocked.Increment(ref progress);
                    Interlocked.Add(ref subcount, sc * autClassSizes[sym]);
                    done[sym, y, x] = true;
                };
            });
            Progress(progress, total, sw.Elapsed + highestTime, subcount);
            Enumerable.Range(0, 4).SelectMany(sym => Enumerable.Range(0, 2 * n).SelectMany(y => Enumerable.Range(0, 2 * n).Select(x => (sym, y, x)))).Where(e => !done[e.sym, e.y, e.x]).OrderBy(e => Guid.NewGuid() /* shuffle */).AsParallel().ForAll(e =>
            {
                int sym = e.sym, y = e.y, x = e.x;
                int i = 1 - n + y, j = 1 - n + Math.Abs(i) + x;
                string fileName = $"{n}-{sym}_{i}_{j}.txt";
                long sc = CountSymmetricPolycubes(matrixReps[sym], new int[3] { i * affine1[sym][0] + j * affine2[sym][0] + biases[sym][0], i * affine1[sym][1] + j * affine2[sym][1] + biases[sym][1], i * affine1[sym][2] + j * affine2[sym][2] + biases[sym][2] });
                done[sym, y, x] = true;
                File.WriteAllText(fileName, $"{sc} {sw.Elapsed + highestTime}");
                Interlocked.Add(ref subcounts[sym], sc);
                Interlocked.Increment(ref progress);
                Interlocked.Add(ref subcount, sc * autClassSizes[sym]);
                Progress(progress, total, sw.Elapsed + highestTime, subcount);
            });
            Console.WriteLine();
            for (int sym = 0; sym < 4; sym++)
            {
                long sc = subcounts[sym];
                Console.WriteLine($"{sc:N0} polycubes fixed under each {descriptions[sym]} * {autClassSizes[sym]} = {autClassSizes[sym] * sc:N0}");
                count += autClassSizes[sym] * sc;
            }
            sw.Stop();
            Console.WriteLine($"{count:N0} total count for non-trivial symmetries for polycubes with {n} cells, {count / 24:N0} free polycubes");
            totalTime += sw.Elapsed + highestTime;
        }
        Console.WriteLine();
        {
            Console.WriteLine("Phase 2/2: Computing trivial symmetries");
            Stopwatch sw = Stopwatch.StartNew();
            int tasks = 4 * (n - filterDepth) - 1;//this is the maximum left stack length, which is the value being filtered to separate work
            long[] subcounts = new long[tasks];
            int progress = 0;
            TimeSpan highestTime = TimeSpan.Zero;
            Enumerable.Range(0, tasks).AsParallel().ForAll(j =>
            {
                string fileName = $"{n}-{filterDepth}-{j}.txt";
                if (File.Exists(fileName))
                {
                    string[] fileParts = File.ReadAllText(fileName).Split(' ');
                    subcounts[j] = long.Parse(fileParts[0]);
                    if (fileParts.Length > 1 && TimeSpan.Parse(fileParts[1]) > highestTime) highestTime = TimeSpan.Parse(fileParts[1]);
                    Interlocked.Add(ref progress, 1);
                }
                else subcounts[j] = -1;
            });
            Progress(progress, tasks, sw.Elapsed + highestTime, count + subcounts.Where(c => c > 0).Sum());
            Enumerable.Range(0, tasks).Where(j => subcounts[j] == -1).OrderBy(j => Guid.NewGuid() /* shuffle */).AsParallel().ForAll(j =>
            {
                string fileName = $"{n}-{filterDepth}-{j}.txt";
                subcounts[j] = CountExtensionsSubset(j, filterDepth);
                File.WriteAllText(fileName, $"{subcounts[j]} {sw.Elapsed + highestTime}");
                Interlocked.Add(ref progress, 1);
                Progress(progress, tasks, sw.Elapsed + highestTime, count + subcounts.Where(c => c > 0).Sum());
            });
            Console.WriteLine();
            long subcount = subcounts.Sum();
            count += subcount;
            sw.Stop();
            Console.WriteLine($"{subcount:N0} polycubes with {n} cells (number of polycubes fixed by trivial symmetry), {subcount / 24:N0} free polycubes");
            totalTime += sw.Elapsed + highestTime;
        }
        Console.WriteLine();
        count /= 24;
        Console.WriteLine($"{count:N0} free polycubes with {n} cells");
        Console.WriteLine($"Total time: {totalTime}");
        return 0;
    }
    private static void Progress(int progress, int total, TimeSpan timeSoFar, long count)
    {
        Console.Write($"\r{count / 24:N0} {timeSoFar} {progress}/{total} {Estimate(progress, total, timeSoFar)}   ");
    }
    private static string Estimate(int progress, int total, TimeSpan time)
    {
        if (progress == 0) return "";
        if (progress == total) return "Completed!                  ";
        double rate = progress / time.TotalSeconds, remaining = (total - progress) / rate;
        return $"ETA {TimeSpan.FromSeconds(remaining)}";
    }
    private static long CountSymmetricPolycubes(int[] linearMap, int[] affineShift)
    {
        byte[,,] adjacencyCounts = new byte[n + 2, 2 * n + 1, 2 * n + 1];
        for (int y = 0; y < 2 * n + 1; y++)
        {
            for (int x = 0; x < 2 * n + 1; x++)
            {
                adjacencyCounts[0, y, x] = 1;
            }
        }
        for (int y = 0; y < n; y++)
        {
            for (int x = 0; x < 2 * n + 1; x++)
            {
                adjacencyCounts[1, y, x] = 1;
            }
        }
        for (int x = 0; x <= n; x++)
        {
            adjacencyCounts[1, n, x] = 1;
        }

        HashSet<(int, int, int)> requiredCells = new();
        Stack<(int, int, int)> extensionStack = new(), recoveryStack = new();
        extensionStack.Push((n, n, 1));
        return CountExtensions(n);

        long CountExtensions(int cellsToAdd)
        {
            cellsToAdd--;
            long count = 0;
            int originalLength = extensionStack.Count;
            while (extensionStack.Count > 0)
            {
                int x, y, z;
                recoveryStack.Push((x, y, z) = extensionStack.Pop());

                bool existingRequirement = requiredCells.Remove((x, y, z));
                if (!existingRequirement)
                {
                    if (cellsToAdd < requiredCells.Count) continue;//number of required cells will only grow, so if already impossible, go to next
                    for (int tempX = x, tempY = y, tempZ = z; ;)//works for general transformations of finite order
                    {
                        (tempX, tempY, tempZ) = (linearMap[0] * tempX + linearMap[1] * tempY + linearMap[2] * tempZ + affineShift[0],
                                                 linearMap[3] * tempX + linearMap[4] * tempY + linearMap[5] * tempZ + affineShift[1],
                                                 linearMap[6] * tempX + linearMap[7] * tempY + linearMap[8] * tempZ + affineShift[2]);
                        if (x == tempX && y == tempY && z == tempZ) break;
                        requiredCells.Add((tempX, tempY, tempZ));
                    }
                }

                if (cellsToAdd >= requiredCells.Count)//if there are too many required cells, then no valid polycubes are possible
                {
                    if (cellsToAdd == 0) count++;
                    else
                    {
                        int innerOriginalLength = extensionStack.Count;
                        if (adjacencyCounts[z, y, x - 1]++ == 0) extensionStack.Push((x - 1, y, z));
                        if (adjacencyCounts[z, y - 1, x]++ == 0) extensionStack.Push((x, y - 1, z));
                        if (adjacencyCounts[z - 1, y, x]++ == 0) extensionStack.Push((x, y, z - 1));
                        if (adjacencyCounts[z, y, x + 1]++ == 0) extensionStack.Push((x + 1, y, z));
                        if (adjacencyCounts[z, y + 1, x]++ == 0) extensionStack.Push((x, y + 1, z));
                        if (adjacencyCounts[z + 1, y, x]++ == 0) extensionStack.Push((x, y, z + 1));

                        count += CountExtensions(cellsToAdd);

                        --adjacencyCounts[z, y, x - 1];
                        --adjacencyCounts[z, y - 1, x];
                        --adjacencyCounts[z - 1, y, x];
                        --adjacencyCounts[z, y, x + 1];
                        --adjacencyCounts[z, y + 1, x];
                        --adjacencyCounts[z + 1, y, x];

                        while (extensionStack.Count != innerOriginalLength)
                            extensionStack.Pop();//should replace this with custom stack to avoid this unnecessary loop
                    }
                }

                if (existingRequirement)
                {
                    requiredCells.Add((x, y, z));
                    break;//this required cell will no longer be available in the extension stack, so no more valid polycubes are possible in this branch
                }

                for (int tempX = x, tempY = y, tempZ = z; ;)
                {
                    (tempX, tempY, tempZ) = (linearMap[0] * tempX + linearMap[1] * tempY + linearMap[2] * tempZ + affineShift[0],
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
    private static long CountExtensionsSubset(int filter, int filterDepth)
    {
        unsafe
        {
            int X = (n + 5) / 4 * ((n + 5) / 4 * 3 - 2);//set X<Y<Z such that aX+bY+cZ = 0 implies a = b = c = 0 or |a|+|b|+|c| > n
            int Y = X + 1;//trivial choice is X = 1, Y = n, Z = n * n. A simple reduction is X = 1, Y = n, Z = n * (n / 2) + (n + 1) / 2
            int Z = X + (n + 5) / 4 * 3;//minimising Z is memory efficient. Unclear if this noticably affects performance. Z ~ 3/16 n^2 is the best I can find for arbitrary n
            byte* byteBoard = stackalloc byte[(n + 2) * Z];//could use ints or shorts as offsets to save memory, but it's faster to directly store the pointers to avoid adding pointer offsets at every lookup
            byte** refStack = stackalloc byte*[(n - 2) * 4];//total length of the two stacks is at most 4n-9. One stack grows from the left, the other stack grows from the right
            *refStack = byteBoard += Z;//seeded with first index of the byte board as the only allowed extension
            for (byte* i = byteBoard + (n + 1) * Z; --i != byteBoard;)
                *i = 255;//the first Z + 1 bytes are disallowed extensions; first Z are less than the minimum, last 1 due to edge case of initial polycube having no neighbours

            long count = CountExtensions(n, refStack + 1, refStack + (n - 2) * 4);
            return count;

            long CountExtensions(int depth, byte** stackTop1, byte** stackTop2)
            {
                long count = 0;
                byte** stackTopOriginal = stackTop1;
                while (stackTop1 != refStack)
                {
                    byte* index = *--stackTop1;
                    byte** stackTopInner = stackTop1;
                    byte* iXm = index - X; if (++*iXm == 0) *stackTopInner++ = iXm;
                    byte* iYm = index - Y; if (++*iYm == 0) *stackTopInner++ = iYm;
                    byte* iZm = index - Z; if (++*iZm == 0) *stackTopInner++ = iZm;
                    byte* iXp = index + X; if (++*iXp == 0) *stackTopInner++ = iXp;
                    byte* iYp = index + Y; if (++*iYp == 0) *stackTopInner++ = iYp;
                    byte* iZp = index + Z; if (++*iZp == 0) *stackTopInner++ = iZp;

                    if (depth == 4) count += CountFinalExtensions(stackTopInner);
                    else if (depth != filterDepth || stackTop1 - refStack == filter) count += CountExtensions(depth - 1, stackTopInner, stackTop2);//if multithreading is not wanted, remove "if (condition)" from this else statement

                    --*iXm;
                    --*iYm;
                    --*iZm;
                    --*iXp;
                    --*iYp;
                    --*iZp;
                    *--stackTop2 = index;//doing this push before the recursion would add one extra unnecessary element to the stack at each level of recursion
                }
                while (stackTop1 != stackTopOriginal)
                    *stackTop1++ = *stackTop2++;
                return count;
            }

            int CountFinalExtensions(byte** stackTop)
            {
                int X2 = X + X, Y2 = Y + Y, Z2 = Z + Z, sYX = Y + X, sZX = Z + X, sZY = Z + Y, dYX = Y - X, dZX = Z - X, dZY = Z - Y;
                int length = (int)(stackTop - refStack), count = length * (length - 1) * (length - 2) / 6, neighbourAdjust = (length << 1) - 511;
                byte** stackTopTemp = stackTop;
                while (stackTopTemp != refStack)
                {
                    byte* i = *--stackTopTemp, p;
                    int neighbours = 0, subcount = 128;
                    p = i - X; if (*p > 127) { count += --*p; neighbours++; subcount += *(i - X2) + *(i - sYX) + *(i - sZX) + *(i + dYX) + *(i + dZX); }
                    p = i - Y; if (*p > 127) { count += --*p; neighbours++; subcount += *(i - Y2) + *(i - sYX) + *(i - sZY) + *(i - dYX) + *(i + dZY); }
                    p = i - Z; if (*p > 127) { count += --*p; neighbours++; subcount += *(i - Z2) + *(i - sZX) + *(i - sZY) + *(i - dZX) + *(i - dZY); }
                    p = i + X; if (*p > 127) { count += --*p; neighbours++; subcount += *(i + X2) + *(i + sYX) + *(i + sZX) + *(i - dYX) + *(i - dZX); }
                    p = i + Y; if (*p > 127) { count += --*p; neighbours++; subcount += *(i + Y2) + *(i + sYX) + *(i + sZY) + *(i + dYX) + *(i - dZY); }
                    p = i + Z; if (*p > 127) { count += --*p; neighbours++; subcount += *(i + Z2) + *(i + sZX) + *(i + sZY) + *(i + dZX) + *(i + dZY); }
                    count += (subcount >> 8) + (neighbours * (neighbours + neighbourAdjust) >> 1);
                }
                while (stackTop != refStack)
                {
                    byte* i = *--stackTop, p;
                    p = i - X; *p |= (byte)(*p >> 4);
                    p = i - Y; *p |= (byte)(*p >> 4);
                    p = i - Z; *p |= (byte)(*p >> 4);
                    p = i + X; *p |= (byte)(*p >> 4);
                    p = i + Y; *p |= (byte)(*p >> 4);
                    p = i + Z; *p |= (byte)(*p >> 4);
                }
                return count;
            }
        }
    }
}

//written by Stanley Dodds (that's me!)
//I can probably supply a full proof of correctness if needed, but I have nothing written down right now - so currently, source: dude, just trust me (or try reading my code)
//Of course, you can compare up to n = 19 with previous results and see that it's definitely doing something right