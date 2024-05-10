# Finite Element Method Implementation for Elastic Bar

## Introduction
This repository contains a generic finite element code developed to solve the governing differential equation for an elastic bar under traction load and constrained at the ends. The code is capable of handling various boundary conditions, variable material properties, and different interpolation methods.

## Features
1. **Boundary Conditions/End Constraints**: Both ends can be constrained by specifying primary variable (Dirichlet/Displacement/Essential), secondary variable or force (Force/Neumann/Natural), and springs (Mixed/Robin).
2. **Variable Material Properties**: The variables T(x), c(x), and AE(x) can vary from a constant to a quadratic function.
3. **Discretization**: The length L and the number of elements can be input values. The domain is discretized into a given number of elements with equal lengths.
4. **Concentrated Load**: Provision to put at least one concentrated load at any given location (excluding the ends).
5. **Interpolation Methods**: Use of Lagrange interpolation or hierarchic shape functions up to quartic order.
6. **Postprocessing**: Represent primary, secondary, and other variables over the domain either continuously or discretely as required.

## Usage
To use the code, follow these steps:
1. Clone the repository: `git clone https://github.com/ppsingh21/FiniteElementMethods.git`
2. Navigate to the project directory: `cd 1hp_code`
3. Run the main code file with appropriate inputs.

## Patch Test
1. Perform patch tests for various cases as described in the problem statement.
2. Compare the solutions obtained with exact solutions.
3. Plot the error in the solution and discuss the results.

## Strain Energy Analysis
1. Plot the strain energy of the exact and finite element solutions against the number of elements in the mesh for all cases.
2. Discuss the results.

## Results Discussion
1. Discuss the results obtained from patch tests and strain energy analysis.
2. Compare solutions obtained with different interpolation methods.
3. Analyze convergence rates and error behaviors.

## Requirements
- Python 3.x
- NumPy
- Matplotlib (for plotting results)

## Contributors
- [Prabal Pratap Singh](https://github.com/ppsingh21)
- [Akshat Kumar](https://github.com/Akshat-2606)

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
