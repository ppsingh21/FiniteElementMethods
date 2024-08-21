# Finite Element Method Implementation

## Elastic Bar
### Introduction
This repository contains a generic finite element code developed to solve the governing differential equation for an elastic bar under traction load and constrained at the ends. The code is capable of handling various boundary conditions, variable material properties, and different interpolation methods.

### Features
1. **Boundary Conditions/End Constraints**: Both ends can be constrained by specifying primary variable (Dirichlet/Displacement/Essential), secondary variable or force (Force/Neumann/Natural), and springs (Mixed/Robin).
2. **Variable Material Properties**: The variables T(x), c(x), and AE(x) can vary from a constant to a quadratic function.
3. **Discretization**: The length L and the number of elements can be input values. The domain is discretized into a given number of elements with equal lengths.
4. **Concentrated Load**: Provision to put at least one concentrated load at any given location (excluding the ends).
5. **Interpolation Methods**: Use of Lagrange interpolation or hierarchic shape functions up to quartic order.
6. **Postprocessing**: Represent primary, secondary, and other variables over the domain either continuously or discretely as required.

### Usage
To use the code, follow these steps:
1. Clone the repository:
   ```bash
   git clone https://github.com/ppsingh21/FiniteElementMethods.git
3. Navigate to the project directory:
   ```bash
   cd 1hp_code
5. Run the main code file with appropriate inputs.

### Patch Test
1. Perform patch tests for various cases as described in the problem statement.
2. Compare the solutions obtained with exact solutions.
3. Plot the error in the solution and discuss the results.

### Strain Energy Analysis
1. Plot the strain energy of the exact and finite element solutions against the number of elements in the mesh for all cases.
2. Discuss the results.

### Results Discussion
1. Discuss the results obtained from patch tests and strain energy analysis.
2. Compare solutions obtained with different interpolation methods.
3. Analyze convergence rates and error behaviors.

### Requirements
- Python 3.x
- NumPy
- Matplotlib (for plotting results)

### Contributors
- [Prabal Pratap Singh](https://github.com/ppsingh21)
- [Akshat Kumar](https://github.com/Akshat-2606)

## One Dimensional Beam Bending
### Introduction
This repository contains a one-dimensional finite element code implemented for solving the beam bending problem using Hermite cubic shape functions. The code is capable of handling various types of transverse loads and applying different combinations of boundary conditions at either end of the beam.

### Problem Details
1. **Uniform Cross Section**: 1 cm x 1 cm
2. **Length of the Beam**: 10 cm
3. **Material Properties**: E = 200 GPa
4. **Transverse Loads**: 
    - Concentrated/point load
    - Uniformly distributed load
    - Point moments at the center of the beam length only
5. **Boundary Conditions**:
    - Specified transverse displacement
    - Specified slope of the transverse displacement
    - Shear force
    - Bending moment

### Finite Element Analysis
Perform the following finite element analyses using the code for 1, 4, 10, 50, and 100 elements:
1. Continuous variation of transverse displacement and its slope.
2. Continuous variation of shear force and bending moment.
3. Bending stress on the topmost line of the beam along its entire length.

### Results Discussion
Discuss the results obtained from the finite element analyses and verify them using Euler-Bernoulli beam theory closed-form solutions.

### Usage
To use the code, follow these steps:
1. Clone the repository:
   ```bash
   git clone https://github.com/ppsingh21/FiniteElementMethods.git
3. Navigate to the project directory:
   ```bash
   git clone https://github.com/ppsingh21/FiniteElementMethods.git
5. Run the main code file with appropriate inputs.

### Requirements
- Python 3.x
- NumPy

### Contributors
- [Prabal Pratap Singh](https://github.com/ppsingh21)
- [Akshat Kumar](https://github.com/Akshat-2606)

### License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
