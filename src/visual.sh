#!/bin/bash
set -e

read -p "Input filename: " base

# compile C++ code
make

# generate images for both individual and accumulate files
for i in $(seq 1 16); do
    echo "Processing ${base}_individual${i}..."
    python3 ../visual/visualizer.py "${base}_individual${i}"

    echo "Processing ${base}_accumulate${i}..."
    python3 ../visual/visualizer.py "${base}_accumulate${i}"
done

echo "All done. Check ../visual/ for images."
