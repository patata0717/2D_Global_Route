# 2D_Global_Route

This is the project 4 for CSR5302, NTHU.

Group 7, Laura, Patata.

## How to run

2D_Global_Route/src$> make

2D_Global_Route/src$> ../bin/PJ4 <input_file> <output_file>

2D_Global_Route/src$> ../bin/PJ4 ../testcase/sample.txt ../output/sample.out

## How to visualize

For all imgae:

16 images to demonstrate accumulation

16 images to display individual nets

2D_Global_Route/src$> ./visual.sh

Input filename? > <output_header>

Input filename? > sample

## Feature

- Sequential method
- Stage 1: Find best 1-steiner by Kahng and Robin method for each net O(N3), which give the smallest wirelength.
- Stage 2: Try to fix overflow by Rip-up Reroute, reroute by Lee algorithm O(MN).

## Experiment Result(for sample)

- Best steiner route

Wirelength = 83

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

![imgur](https://ibb.co/23DY6zFj)

![imgur](https://ibb.co/9k6Nwcsw)
