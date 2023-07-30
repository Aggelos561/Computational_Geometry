# Polygonization and Optimization Algorithms

## Introduction

This project aims to find the optimal polygonization of point sets using various algorithms. It provides implementations of incremental, convex hull, and subdivision polygonization algorithms, along with local search and simulated annealing optimization algorithms.

## Compilation

To compile the project, follow these steps in the `src/` directory:

1. Generate the CMakeLists.txt for the 'evaluate' target:
   ```
   cgal_create_CMakeLists -s evaluate
   ```

2. Configure the build environment:
   ```
   cmake -DCMAKE_BUILD_TYPE=Release .
   ```

3. Build the project:
   ```
   make
   ```

## Algorithm Execution

To execute a specific algorithm and optimization, use the following command:

```bash
$ ./evaluate -i (point set path) -o (output file) -algorithm (optimization_algorithm) -L (L_value) -max|-min -annealing (annealing_method) -algorithm_initial (initial_algorithm) -initialization (initialization_option) -m (subdivisions) -edge_selection (edge_selection_method)
```

### Command-line Parameters

- `(point set path)`: `.instance` file to be processed.
- `(output file)`: File to record the data and results
- `optimization_algorithm`: The optimization algorithm to use ('local_search' or 'simulated_annealing').
- `L_value`: The value of the parameter 'L' (depends on the algorithm).
- `-max` or `-min`: Choose between maximizing or minimizing the objective function.
- `annealing_method`: The annealing method for simulated annealing ('local', 'global', or 'subdivision').
- `initial_algorithm`: The initial polygonization algorithm ('incremental' or 'convex_hull').
- `initialization_option`: Applicable only for the incremental algorithm. Choose between '1a', '1b', '2a', or '2b'. Subdivision is only allowed with '1a' and '1b'.
- `subdivisions`: The number of subdivisions (only relevant for subdivision simulated annealing).
- `edge_selection_method`: The method for selecting initial edges (1 for random, 2 for minimum, 3 for maximum). If not provided, the optimization algorithm decides.

### Example Execution

For example, to run the optimization algorithm 'simulated_annealing' with the following parameters:

- Input file: `uniform-0000100-1.instance`
- Output file: `output.txt`
- 'L' value: 10000
- Maximize the objective function
- Annealing method: 'global'
- Initial polygonization algorithm: 'incremental' with initialization '1a'
- Use random edge selection

The command will be:

```bash
$ ./evaluate -i ../instances_test/uniform-0000100-1.instance -o output.txt -algorithm simulated_annealing -L 10000 -max -annealing global -algorithm_initial incremental -initialization 1a -edge_selection 1
```

## Scoring with Preprocessing

To run the scoring tables with preprocessing enabled, use the following command:

```bash
$ ./evaluate -i (point set path) -o (output file) -preprocess -scores
```

### Command-line Parameters

- `(point set path)`: Directory containing all the `.instance` files to be processed.
- `(output file)`: File to record the data and results (result table).
- `-preprocess`: This parameter is optional. If executed with `-preprocess`, heuristic algorithms will run on each file of the input directory, and an average for the algorithm parameter 'L' will be determined based on the point set size. For larger point sets (greater than 1000 points), the algorithm creates 2 subsets with random points and performs tests for 'L' to find an optimal value in a relatively short period of time.

### Visualization

The project includes a Python script `graph.py` for visualizing the resulting polygons. 

## Results

The project generates two result files:

1. `resultsWithoutPreprocess.txt`: Contains the scoring tables without any preprocessing.
2. `resultsWithPreprocess.txt`: Contains the scoring tables after preprocessing if enabled using the `-preprocess` flag.

## Notes on Algorithm Performance

1. In most cases, the Subdivision algorithm yields better results for the minimum score compared to other algorithms.

2. For the maximum score, Spatial Subdivision and Convex+Local algorithms provide better results, with Convex+Local showing improvement as the number of points increases but failing to meet the cut-off time for inputs greater than or equal to 10,000 points.

3. Subdivision exhibits the shortest execution time, while Convex+Local has the longest execution time.

4. Incr+Global+Local and Incr+Local algorithms fall in between regarding execution time, with Incr+Global+Local outperforming in final results due to the addition of Global Transitions, which enhances the final area for almost every input.

5. All selected algorithms, except Convex+Local, manage to produce results within the cut-off time for inputs up to 100,000 points.

6. Running algorithms with preprocessing generally yields better results, particularly in minimization cases, but also in many cases of maximization.

7. The local transitions of the simulated annealing, due to the small local changes and the use of a kd-tree for intersection checks, allow for larger values of 'L' compared to the global transitions.

8. The execution time of Spatial Subdivision increases linearly due to dividing the point set into subsets, resulting in the algorithm being executed on the same number of points, but with an increasing number of subsets.

9. The simulated annealing algorithm provides better approximations for the area, while being satisfactory in terms of execution time, especially when combined with the Subdivision annealing. Local transitions yield faster results but less optimized, whereas global transitions are slower but produce more optimal results. The Spatial Subdivision combines global transitions on small polygons and local transitions on the final polygon, resulting in a satisfactory outcome in both time and approximation.
