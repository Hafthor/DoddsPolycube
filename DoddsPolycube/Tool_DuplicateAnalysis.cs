namespace DoddsPolycube;

using System.Runtime.InteropServices;
using Num = ulong;

// Cache analysis: measure duplicate board states at each depth
// =============================================================

public partial class Program {

    /// <summary>
    /// Measures duplicate rates using two hash strategies:
    /// 1. Full board hash (exact state)
    /// 2. Extension-set hash (just the positions on the stack, sorted)
    /// Strategy 2 asks: "do different placements of cells ever produce
    /// the same set of extension positions?" If so, we could cache.
    /// Single-threaded, all filters, depths 4-8.
    /// </summary>
    private static unsafe void MeasureDuplicates() {
        Console.WriteLine($"Duplicate analysis: N={N}, FilterDepth={FilterDepth}");
        Console.WriteLine();

        int boardSize = (N + 2) * Z;
        int stackSize = (N - 2) * 4;

        byte* byteBoard = stackalloc byte[boardSize];
        byte** refStack = stackalloc byte*[stackSize];
        byte* boardStart = byteBoard + Z;
        *refStack = boardStart;
        for (byte* i = boardStart + (N + 1) * Z; --i != boardStart;)
            *i = 255;

        int trackFromDepth = 4;
        int trackToDepth = Math.Min(FilterDepth, 8);
        long[] totalVisits = new long[trackToDepth + 1];

        // Track duplicates via board hash and extension-set hash
        var boardHashSets = new Dictionary<int, HashSet<ulong>>();
        var extHashSets = new Dictionary<int, HashSet<ulong>>();
        for (int d = trackFromDepth; d <= trackToDepth; d++) {
            boardHashSets[d] = new HashSet<ulong>();
            extHashSets[d] = new HashSet<ulong>();
        }

        var sw = System.Diagnostics.Stopwatch.StartNew();
        Num totalCount = 0;
        int totalFilters = MaxLeftStackLen + 1;
        for (int filter = 0; filter < totalFilters; filter++) {
            Num cnt = CountWithTracking(
                N, refStack + 1, refStack + stackSize, filter,
                boardStart, boardSize, refStack,
                totalVisits, boardHashSets, extHashSets, trackFromDepth, trackToDepth);
            totalCount += cnt;
            Console.Write($"\r  filter {filter + 1}/{totalFilters} ({sw.Elapsed:mm\\:ss})");
        }
        Console.WriteLine($"\r  all {totalFilters} filters done in {sw.Elapsed:mm\\:ss}.          ");
        Console.WriteLine($"Total count: {totalCount:N0}");
        Console.WriteLine();
        Console.WriteLine("Depth  TotalVisits      BoardUnique      ExtSetUnique     BoardHitRate  ExtSetHitRate");
        Console.WriteLine("-----  -----------      -----------      ------------     ------------  -------------");

        for (int d = trackFromDepth; d <= trackToDepth; d++) {
            long total = totalVisits[d];
            long boardUniq = boardHashSets[d].Count;
            long extUniq = extHashSets[d].Count;
            if (total == 0) continue;
            double boardHit = total > boardUniq ? (double)(total - boardUniq) / total : 0;
            double extHit = total > extUniq ? (double)(total - extUniq) / total : 0;
            Console.WriteLine($"  {d,3}  {total,15:N0}  {boardUniq,15:N0}  {extUniq,15:N0}     {boardHit,10:P1}   {extHit,10:P1}");
        }

        Console.WriteLine();
        Console.WriteLine("Depth  ExtSetUnique   EstCacheMemory  EstSavings");
        Console.WriteLine("-----  ------------   --------------  ----------");
        for (int d = trackFromDepth; d <= trackToDepth; d++) {
            long total = totalVisits[d];
            long extUniq = extHashSets[d].Count;
            if (extUniq == 0) continue;
            long memBytes = extUniq * 48;
            string memStr = memBytes < 1024 * 1024 ? $"{memBytes / 1024.0:F0} KB"
                : memBytes < 1024L * 1024 * 1024 ? $"{memBytes / (1024.0 * 1024):F1} MB"
                : $"{memBytes / (1024.0 * 1024 * 1024):F1} GB";
            double savings = total > 0 ? (double)(total - extUniq) / total : 0;
            Console.WriteLine($"  {d,3}  {extUniq,12:N0}   {memStr,14}  {savings,8:P1}");
        }

        static ulong HashBoard(byte* board, int size) {
            ulong hash = 14695981039346656037UL;
            byte* end = board + size;
            while (board + 8 <= end) {
                hash ^= *(ulong*)board;
                hash *= 1099511628211UL;
                board += 8;
            }
            while (board < end) {
                hash ^= *board++;
                hash *= 1099511628211UL;
            }
            return hash;
        }

        // Hash the extension set: sort the stack offsets and hash them.
        // This is independent of placement order — only depends on WHICH
        // positions are available extensions.
        static ulong HashExtensionSet(byte** stackBase, byte** stackTop, byte* boardBase) {
            int len = (int)(stackTop - stackBase);
            // Collect offsets
            Span<int> offsets = stackalloc int[len];
            for (int i = 0; i < len; i++)
                offsets[i] = (int)(stackBase[i] - boardBase);
            offsets.Sort();
            // Hash sorted offsets
            ulong hash = 14695981039346656037UL;
            for (int i = 0; i < len; i++) {
                hash ^= (ulong)offsets[i];
                hash *= 1099511628211UL;
            }
            // Also include depth (number of placed cells) to distinguish
            // same extension sets at different depths
            return hash;
        }

        Num CountWithTracking(
            int depth, byte** stackTop1, byte** stackTop2, int filter,
            byte* boardBase, int bSize, byte** sBase,
            long[] visits,
            Dictionary<int, HashSet<ulong>> bSets,
            Dictionary<int, HashSet<ulong>> eSets,
            int tFrom, int tTo) {

            Num cnt = 0;
            byte** stackTopOriginal = stackTop1;
            while (stackTop1 != sBase) {
                byte* index = *--stackTop1;
                if (depth != FilterDepth || stackTop1 - sBase == filter) {
                    byte** stackTopInner = stackTop1;

                    if (++*(index - Z) == 0) *stackTopInner++ = index - Z;
                    if (++*(index - Y) == 0) *stackTopInner++ = index - Y;
                    if (++*(index - X) == 0) *stackTopInner++ = index - X;
                    if (++*(index + X) == 0) *stackTopInner++ = index + X;
                    if (++*(index + Y) == 0) *stackTopInner++ = index + Y;
                    if (++*(index + Z) == 0) *stackTopInner++ = index + Z;

                    if (depth >= tFrom && depth <= tTo) {
                        visits[depth]++;
                        bSets[depth].Add(HashBoard(boardBase, bSize));
                        eSets[depth].Add(HashExtensionSet(sBase, stackTopInner, boardBase));
                    }

                    if (depth == 4) {
                        byte** stackTop = stackTopInner, stackTopTemp = stackTopInner;
                        int length = (int)(stackTop - sBase), lengthPlus = (length << 1) - 511;
                        cnt += (Num)(length * (length - 1) * (length - 2) / 6);
                        for (; stackTopTemp != sBase;) {
                            byte* i = *--stackTopTemp;
                            int neighbours = 0, subCount = 128, localCount = 0;
                            if (*(i - Z) > 127) { localCount += --*(i - Z); subCount += *(i-Z-Z)+*(i-Z-X)+*(i-Z-Y)+*(i-Z+X)+*(i-Z+Y); neighbours++; }
                            if (*(i - Y) > 127) { localCount += --*(i - Y); subCount += *(i-Y-Y)+*(i-Y-X)+*(i-Y-Z)+*(i-Y+X)+*(i-Y+Z); neighbours++; }
                            if (*(i - X) > 127) { localCount += --*(i - X); subCount += *(i-X-X)+*(i-X-Y)+*(i-X-Z)+*(i-X+Y)+*(i-X+Z); neighbours++; }
                            if (*(i + X) > 127) { localCount += --*(i + X); subCount += *(i+X+X)+*(i+X+Y)+*(i+X+Z)+*(i+X-Y)+*(i+X-Z); neighbours++; }
                            if (*(i + Y) > 127) { localCount += --*(i + Y); subCount += *(i+Y+Y)+*(i+Y+X)+*(i+Y+Z)+*(i+Y-X)+*(i+Y-Z); neighbours++; }
                            if (*(i + Z) > 127) { localCount += --*(i + Z); subCount += *(i+Z+Z)+*(i+Z+X)+*(i+Z+Y)+*(i+Z-X)+*(i+Z-Y); neighbours++; }
                            cnt += (Num)localCount + (Num)((subCount >> 8) + (neighbours * (neighbours + lengthPlus) >> 1));
                        }
                        while (stackTop != sBase) {
                            byte* i = *--stackTop;
                            *(i-Z)|=(byte)(*(i-Z)>>4); *(i-Y)|=(byte)(*(i-Y)>>4); *(i-X)|=(byte)(*(i-X)>>4);
                            *(i+X)|=(byte)(*(i+X)>>4); *(i+Y)|=(byte)(*(i+Y)>>4); *(i+Z)|=(byte)(*(i+Z)>>4);
                        }
                    } else {
                        cnt += CountWithTracking(
                            depth - 1, stackTopInner, stackTop2, filter,
                            boardBase, bSize, sBase,
                            visits, bSets, eSets, tFrom, tTo);
                    }

                    --*(index - Z); --*(index - Y); --*(index - X);
                    --*(index + X); --*(index + Y); --*(index + Z);
                }
                *--stackTop2 = index;
            }
            while (stackTop1 != stackTopOriginal)
                *stackTop1++ = *stackTop2++;
            return cnt;
        }
    }
}

