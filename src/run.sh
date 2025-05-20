#!/bin/bash

# Compile the C++ code
make

# Run
for i in $(seq 1 16); do
    # Run the executable with input, output, and net index
    ../bin/PJ4 ../testcase/sample.txt ../output/sample${i}.out 1 $i

    # Generate visualization
    python3 ../visual/visualizer.py sample${i}
done

echo "All done. Check ../visual/ for images."
