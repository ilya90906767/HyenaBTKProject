#!/bin/bash

input_dir=/mnt/data/ibetyaev/hyena/window_filtered/
output_filtered_dir=/mnt/data/ibetyaev/hyena/window_filtered/trees/

# Create the output directory if it doesn't exist
mkdir -p "$output_filtered_dir"

# Loop through all files in the input directory
for file in "$input_dir"*.phy; do
  # Get the filename without the directory and extension
  filename="${file##*/}"
  filename="${filename%.phy}"

  output_subdir="$output_filtered_dir/${filename}_raxml"
  mkdir -p "$output_subdir"

  # Run RAxML-NG with BINGAMMA model and outgroup
  raxml-ng --search --model BINGAMMA --threads 1 --msa "$file" --prefix "$(dirname "$output_subdir")"
done