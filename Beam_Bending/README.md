# One Dimensional Finite Element Code for Beam Bending Problem

## Introduction
This repository contains a one-dimensional finite element code implemented for solving the beam bending problem using Hermite cubic shape functions. The code is capable of handling various types of transverse loads and applying different combinations of boundary conditions at either end of the beam.

## Problem Details
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

## Finite Element Analysis
Perform the following finite element analyses using the code for 1, 4, 10, 50, and 100 elements:
1. Continuous variation of transverse displacement and its slope.
2. Continuous variation of shear force and bending moment.
3. Bending stress on the topmost line of the beam along its entire length.

## Results Discussion
Discuss the results obtained from the finite element analyses and verify them using Euler-Bernoulli beam theory closed-form solutions.

## Usage
To use the code, follow these steps:
1. Clone the repository: `git clone https://github.com/ppsingh21/FiniteElementMethods.git`
2. Navigate to the project directory: `cd Beam_Bending`
3. Run the main code file with appropriate inputs.

## Requirements
- Python 3.x
- NumPy

## Contributors
- [Prabal Pratap Singh](https://github.com/ppsingh21)
- [Akshat Kumar](https://github.com/Akshat-2606)

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
