#!/bin/bash

# Directory containing the files
DIR="/mnt/data/ibetyaev/hyena/aligned"
DIR_count="/mnt/data/ibetyaev/hyena/temp/consensus/SRR9914655.fastp_merged_only_GCA_003009895_bwa_mem_sorted_q24_RG_MD_consensus/"
to_move="/mnt/data/ibetyaev/hyena/window"

# Initialize count of files with exactly 6 lines containing '>'
count=0

# Iterate through all files in the directory
for file in "$DIR"/*
do
    if [ -f "$file" ]; then
        # Count the number of lines containing '>'
        line_count=$(grep -c '>' "$file")
        
        # Check if the count is exactly 6
        if [ "$line_count" -eq 8 ]; then
            count=$((count + 1))
            cp "$file" "$to_move"
        fi
    fi
done

total_files_in_count_dir=0

for file in "$DIR_count"/* 
do 
    if [ -f "$file" ]; then
        total_files_in_count_dir=$((total_files_in_count_dir + 1))
    fi
done
# Output the result
echo "Number of files with exactly 9 lines containing '>': $count"
echo "Percente is $count/$total_files_in_count_dir"