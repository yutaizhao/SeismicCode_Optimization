#!/bin/bash

# Create a temporary directory to store intermediate files
tmpdir=$(mktemp -d)

# Loop 1000 times
for ((i=1; i<=101; i++))
do
    # Run the program and calculate the mean of the 2nd column
    output=$(./top-stencil | awk '{ sum += $2; n++ } END { if (n > 0) print sum / n }')
    
    # Save the mean to a file
    echo "$output" >> "$tmpdir/means.txt"
done

# Calculate mean, median, min, max, and stddev of the means
mean=$(awk '{ sum += $1; n++ } END { if (n > 0) print sum / n }' "$tmpdir/means.txt")
median=$(sort -n "$tmpdir/means.txt" | awk 'NR == FNR { a[NR] = $1; next } END { if (NR % 2 == 0) print (a[NR / 2] + a[NR / 2 + 1]) / 2; else print a[(NR + 1) / 2] }')
min=$(awk 'NR == 1 { min = $1 } $1 < min { min = $1 } END { print min }' "$tmpdir/means.txt")
max=$(awk 'NR == 1 { max = $1 } $1 > max { max = $1 } END { print max }' "$tmpdir/means.txt")
stddev=$(awk -v mean="$mean" '{ sum += ($1 - mean)^2; n++ } END { if (n > 1) print sqrt(sum / (n - 1)) }' "$tmpdir/means.txt")

# Print the results
echo "Minimum of means: $min"
echo "Maximum of means: $max"
echo "Mean of means: $mean"
echo "Median of means: $median"
echo "Standard deviation of means: $stddev"

# Clean up
rm -rf "$tmpdir"

