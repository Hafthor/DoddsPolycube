# Plan for Computing Free Polycubes at N=23

## Current Configuration
- **N=23**, **FilterDepth=8**, **Num=ulong** (sufficient for N≤23)
- **MaxLeftStackLen=58** → **59 filter values**
- **Mode**: `--fastcheckpoint` (Exp8) — checkpointable + fast unsafe inner path
- **Threads**: 8 (one per M2 Max performance core)

## How FD=8 Was Determined

### Method 1: Full benchmark at N=16 (`--fdbenchmark -8`)
Tested all FD values 5–14 with the complete computation:
```
FD=5: 62.5s   FD=6: 54.2s ◄   FD=7: 58.6s   FD=8: 54.6s ◄
FD=9: 54.7s   FD=10: 55.2s    FD=11: 57.8s   FD=12: 57.3s
```
**Plateau from FD=6 to FD=10**, all ~13% faster than FD=5.

### Method 2: Quick benchmark at N=16 (`--fdquick -8`)
Single-threaded filter-0-only test (1 minute total, isolates inner path speed):
```
FD=5: 9.2s   FD=6: 3.7s   FD=7: 2.7s ◄   FD=8: 2.7s ◄
FD=9: 3.5s   FD=10: 4.0s  FD=11: 4.7s     FD=12: 5.9s
```
**FD=7-8 optimal for inner computation speed.**

### Why FD=8 for N=23
- **Inner speed**: FD=7 and FD=8 are tied for fastest (quick benchmark)
- **Load balance**: FD=8 → 59 filters / 8 threads = 7.375/thread (excellent)
- **Checkpoint granularity**: 59 trivial result files — fine granularity for resume
- **Outer depth**: 14 managed levels — plenty of checkpoint points
- **Inner depth**: 3 recursion levels + Work4 — optimal unsafe depth

FD=7 (63 filters, 7.875/thread) would also be excellent. The difference is within noise.

## To Determine FD Empirically at Larger N

If you want to validate at a higher N before committing:

1. **Set N=18 or N=19** in Program.cs (takes ~hours instead of months)
2. Run `--fdbenchmark -8` for a full multi-threaded test
3. Or run `--fdquick -8` for a faster single-threaded relative comparison
4. The winning FD at N=18 will be the same as N=23 (plateau is N-independent)

## Production Run Command

```bash
cd /Users/hafthor/Desktop/personal/exp/cs/DoddsPolycube
dotnet run -c Release --project DoddsPolycube -- --fastcheckpoint -8 --checkpointlog
```

This will:
1. **Phase 1**: Compute nontrivial symmetries (~minutes)
2. **Phase 2**: Compute trivial symmetry (the long part)
   - Runs 8 parallel paths, each processing ~7-8 filter values
   - Saves per-filter results to `trivial_23_58_*.txt` as each filter completes
   - Checkpoints to `checkpoint_23_fd8_8p_*.txt` every 60 seconds
   - Appends to `checkpoint_23_fd8.log` on each checkpoint save (proof of work)
   - Displays spinner while running

### Checkpoint Log (`--checkpointlog`)

Adding `--checkpointlog` appends a line to `checkpoint_23_fd8.log` each time any
path saves a checkpoint (~every 60 seconds per path). The log is append-only and
survives restarts, providing a continuous proof-of-work trail. Each line contains:

```
2026-03-04T03:13:53.227Z checkpoint_23_fd8_8p_2.txt depth=9 stackPtr=4 callStackPtr=28 filters=4/59 runningTotal=164534803285
```

| Field | Meaning |
|-------|---------|
| Timestamp | ISO 8601 UTC — when the checkpoint was saved |
| Checkpoint file | Which path wrote this entry |
| depth | Current recursion depth in the outer loop |
| stackPtr | Current stack pointer (lower = deeper in tree) |
| callStackPtr | Call stack depth (tracks nested levels) |
| filters | Completed filters / total for this path |
| runningTotal | Sum of all filter counts so far (monotonically increasing) |
| QUITING | Appended if the save was triggered by Ctrl+C |

Useful for:
- **Proof of work**: Timestamped evidence of computation progress
- **Diagnostics**: If results seem wrong, the log shows the trajectory of running totals
- **Progress estimation**: Compare runningTotal growth rate over time to estimate completion
- **Interruption forensics**: QUITING entries show when/why the process stopped

### Interruption & Resume
- **Ctrl+C**: Saves checkpoint, exits cleanly
- **Restart**: Same command auto-loads checkpoints and completed filters
- **Checkpoint files**: Include N, FD, and path count — safe against recompilation with different settings
- **Trivial files**: Include N and MaxLeftStackLen (which encodes FD) — also safe
- **Log file**: Append-only, accumulates across restarts

### Estimated Runtime
- N=21 took ~X hours. N=23 is roughly **50-100x** longer (branching factor ~7-8 per level, two extra levels).
- With FD=8 (13% faster than FD=5): expect **weeks to months** on M2 Max.

### Monitoring

```bash
# Count completed filters (out of 59)
ls -1 trivial_23_58_*.txt 2>/dev/null | wc -l

# Check checkpoint freshness (are paths still making progress?)
ls -lt checkpoint_23_fd8_8p_*.txt

# Tail the checkpoint log for live progress
tail -f checkpoint_23_fd8.log

# Count total checkpoint log entries (proof of work)
wc -l checkpoint_23_fd8.log
```

