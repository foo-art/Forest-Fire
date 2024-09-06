# Forest Fire

## ```burning_trees.cpp```
### Model Overview
Implementation of the forest fire model in C++ while using OpenMP to parallelise a single run of the model, i.e. it should be possible to run the code for a single random starting grid, grid size N, and
probability p, while using OpenMP. The code output the number of steps before the fire stops burning, whether or not the fire reached the bottom of the grid, and the time taken to run a given simulation, in a format suitable for further analysis. The code are capable of running the forest fire model for M randomly generated starting grids (e.g. using M different random seeds), and averaging over the results. 

### Model Convergence
we talked about convergence with respect to a) the number of random starting grids, M (i.e. the number of repeat runs) and b) the size of the grid, N. For example, as you increase the grid size, the behaviour of the model will change. We expect that at some point these changes will no longer be significant, i.e. for a large enough value of N, further increases will result in only negligible changes. In other words, the model has converged. For this project, we will focus only on convergence with respect to the number of repeat runs, M.

### Model Extension : Wind Direction
In the next part of the project, we will add ‘wind’. This could be done in various ways, but we will use the following rules:
• The wind can only travel up, down, left or right.
• The fire can travel twice as far in the direction of the wind, i.e. the fire can ‘jump’ a cell in the direction of the wind, irrespective of whether or not there is a tree in between.
• The movement of the fire is unaffected in other directions, i.e. continues to travel between nearest neighbours only.

### OpenMP Performance Analysis
The final aspect of the project is to explore the parallel performance of the model using BlueCrystal4. For this part of the project, we are not interested in the behaviour of the model itself, only the time
taken to complete a run, and how this varies with the number of threads. Since the performance is intrinsically linked to the problem size, we study the timing data for three different values of N = 50, 100, 500.

## ```burning_forest_animated.ipynb```
### Visual Simulation
Animation to simulate the dynamic process of forest fire algorithm. I separated the task for the simulation into two parts.
- Part 1 : Class ```forest``` that stores attributes such as matrix to represent the forest and the dimension of the forest. 
- Part 2 : Example for the visual simulation of forest fire algorithm using object ```forest```.

