#!/usr/bin/env bash
# Usage: bench.sh <bam> <outprefix> [jar] [timeout_sec]
# Example: bench.sh low.bam low GeMoSeq-1.2.3.jar 1800
set -u
BAM="$1"
PREFIX="$2"
JAR="${3:-GeMoSeq-1.2.3.jar}"
TIMEOUT="${4:-0}"
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$DIR"
PS1_WIN="$(cygpath -w "$DIR/bench.ps1")"
echo "[bench.sh] start $(date '+%F %T')  bam=$BAM prefix=$PREFIX jar=$JAR timeout=$TIMEOUT"
powershell -NoProfile -ExecutionPolicy Bypass -File "$PS1_WIN" -Bam "$BAM" -Prefix "$PREFIX" -Jar "$JAR" -TimeoutSec "$TIMEOUT"
echo "[bench.sh] end   $(date '+%F %T')  result: $(cat "${PREFIX}.bench.txt" 2>/dev/null)"
