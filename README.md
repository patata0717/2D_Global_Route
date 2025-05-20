# 2D_Global_Route

This is the project 4 for CSR5302, NTHU.

Group 7, Laura, Patata.

## How to run

2D_Global_Route/src$> ../bin/PJ4 <input_file> <output_file> <bottom_layer> <top_layer> 

2D_Global_Route/src$> ../bin/PJ4 ../testcase/sample.txt ../output/sample.out 1 16

<bottom_layer> and <top_layer> is for printing convinience, so set it to 1 and Top_layer, in sample, 16.

## Feature

- Sequential method
- Stage 1: Find best 1-steiner by Kahng and Robin method for each net, which give the smallest wirelength.
- Stage 2: Try to fix overflow by Rip-up Reroute.


## Experiment Result

- Random seed = 12345

Iterations cost = 54, Wirelength = 96

- Random seed = 54321

Iterations cost = 20, Wirelength = 95

- Random seed = 48763

Iterations cost = 110, Wirelength = 96

- True Random

Iterations cost = 171, Wirelength = 97

Iterations cost = 72, Wirelength = 96

Iterations cost = 33, Wirelength = 96
