#!/bin/bash

tmpdir=$(mktemp -d)
for ((i=1; i<=101; i++))
do
    output=$(./top-stencil | awk '{ sum += $2; n++ } END { if (n > 0) print sum / n }')
    
    echo "$output" >> "$tmpdir/means.txt"
done

mean=$(awk '{ sum += $1; n++ } END { if (n > 0) print sum / n }' "$tmpdir/means.txt")
median=$(sort -n "$tmpdir/means.txt" | awk 'NR == FNR { a[NR] = $1; next } END { if (NR % 2 == 0) print (a[NR / 2] + a[NR / 2 + 1]) / 2; else print a[(NR + 1) / 2] }')
min=$(awk 'NR == 1 { min = $1 } $1 < min { min = $1 } END { print min }' "$tmpdir/means.txt")
max=$(awk 'NR == 1 { max = $1 } $1 > max { max = $1 } END { print max }' "$tmpdir/means.txt")
stddev=$(awk -v mean="$mean" '{ sum += ($1 - mean)^2; n++ } END { if (n > 1) print sqrt(sum / (n - 1)) }' "$tmpdir/means.txt")

echo "Minimum of means: $min"
echo "Maximum of means: $max"
echo "Mean of means: $mean"
echo "Median of means: $median"
echo "Standard deviation of means: $stddev"
echo "$min & $max & $mean & $median & $stddev"

rm -rf "$tmpdir"

