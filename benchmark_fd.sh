#!/bin/zsh
# Benchmark different FilterDepth values for the current N
# Usage: ./benchmark_fd.sh [N] [threads]

set -e

N=${1:-16}
THREADS=${2:-8}
PROJECT_DIR="/Users/hafthor/Desktop/personal/exp/cs/DoddsPolycube"
PROGRAM_CS="$PROJECT_DIR/DoddsPolycube/Program.cs"

echo "=== FilterDepth Benchmark for N=$N, $THREADS threads ==="
echo ""

# Save original values
ORIG_N=$(grep "private const int N = " "$PROGRAM_CS" | sed 's/.*N = //;s/;.*//')
ORIG_FD=$(grep "private const int FilterDepth = " "$PROGRAM_CS" | sed 's/.*FilterDepth = //;s/;.*//')

# Set N
sed -i '' "s/private const int N = [0-9]*/private const int N = $N/" "$PROGRAM_CS"

RESULTS=""

# Test FD from 5 to min(N-2, N/2+4)
MAX_FD=$(( N > 14 ? N/2+4 : N-2 ))
for FD in $(seq 5 $MAX_FD); do
    FILTERS=$(( 4 * (N - FD) - 1 ))
    if (( FILTERS < THREADS )); then
        echo "FD=$FD: skipping ($FILTERS filters < $THREADS threads)"
        continue
    fi
    
    sed -i '' "s/private const int FilterDepth = [0-9]*/private const int FilterDepth = $FD/" "$PROGRAM_CS"
    
    if ! dotnet build "$PROJECT_DIR/DoddsPolycube/DoddsPolycube.csproj" -c Release 2>&1 | tail -1 | grep -q "succeeded"; then
        echo "  BUILD FAILED, skipping"
        continue
    fi
    
    echo -n "FD=$FD ($FILTERS filters): "
    OUTPUT=$(dotnet run -c Release --project "$PROJECT_DIR/DoddsPolycube" -- --fdbenchmark -$THREADS 2>&1)
    TIME=$(echo "$OUTPUT" | grep "^→" | sed 's/.*time=//')
    COUNT=$(echo "$OUTPUT" | grep "Total trivial count" | sed 's/.*: //')
    echo "$TIME  (count=$COUNT)"
    RESULTS="$RESULTS\nFD=$FD  $TIME  count=$COUNT"
done

echo ""
echo "=== Summary ==="
echo -e "$RESULTS"

# Restore original values
sed -i '' "s/private const int N = [0-9]*/private const int N = $ORIG_N/" "$PROGRAM_CS"
sed -i '' "s/private const int FilterDepth = [0-9]*/private const int FilterDepth = $ORIG_FD/" "$PROGRAM_CS"
echo ""
echo "Restored N=$ORIG_N, FilterDepth=$ORIG_FD"



