# Decomposition of Symmetrical Classes of Central Configurations

## Authors
**Marcelo P. Santos**  
**Leon D. da Silva**

## Overview
This repository contains the SageMath and Jupyter Notebook implementations for the paper *Decomposition of Symmetrical Classes of Central Configurations*. The research focuses particularly on symmetric configurations of nested polyhedra such as cubes, octahedra, and tetrahedra.

By utilizing the representation theory of finite groups, this project simplifies and decomposes the equations governing central configurations. The code provided here enables symbolic analysis of these configurations, allowing for a deeper exploration of the mass distribution and geometric symmetry required for each configuration.

## Directory Structure

- **Nested Cube/**: Contains scripts for analyzing central configurations of nested cubes.
- **Nested Octahedra/**: Contains scripts for analyzing blocks of nested octahedra.
- **Nested Tetrahedra/**: Contains scripts for analyzing blocks of nested tetrahedra.
- **Example_1/**: Contains scripts for the example 1.
- **util.sage**: Contains utility functions used across different notebooks for symbolic computations in SageMath.

Each directory corresponds to a specific polyhedral configuration and contains:
  - **Jupyter Notebooks**: These `.ipynb` files walk through detailed calculations for each block, following the structure of the paper.
  - **SageMath Scripts**: Supplementary `.sage` scripts are included to support the mathematical computations required in each analysis.

## Requirements

To run the code in this repository, you will need:
- **SageMath 10.4**
- **Jupyter Notebook**

### Installing SageMath 10.4

 Visit the official SageMath page at [SageMath Installation Guide](https://doc.sagemath.org/html/en/installation/index.html) 

## Usage Guide

Each notebook is self-contained and follows the structure of the corresponding sections in the research paper:
- **Nested Cube**: Analysis of central configurations in nested cube structures.
- **Nested Octahedra**: Analysis of central configurations in nested octahedral structures.
- **Nested Tetrahedra**: Analysis of central configurations in nested tetrahedral structures.
- **blocks**/: Contains the matrices resulting from the decomposition in blocks

### Utility Functions
The `util.sage` file provides utility functions used throughout the notebooks for symbolic algebra and matrix manipulations.



Vou refazer a tabela de funções com base no conteúdo correto do arquivo `util.sage` que você forneceu. Aqui está a versão atualizada:


| Function Name                    | Description                                                                   |
| -------------------------------- | ----------------------------------------------------------------------------- |
| `simplify_expression()`          | Simplifies a given polynomial expression by combining like terms iteratively. |
| `count_sign_changes()`           | Counts the number of sign changes in a list of coefficients.                  |
| `sturm()`                        | Uses Sturm's theorem to count the number of real roots of a polynomial in an interval. |
| `expr_to_poly()`                 | Converts a symbolic expression to a polynomial in a specified ring.           |
| `polynomial_to_dict()`           | Converts a polynomial to a dictionary representation with keys as multi-degrees. |
| `moebius_coefficient()`          | Calculates the Möbius coefficient for a given polynomial expression.          |
| `moebius_sign()`                 | Determines the sign of Möbius coefficients for a given polynomial.            |
| `quadratic_extrema_on_interval()`| Finds the minimum and maximum of a quadratic polynomial on a given interval.  |


## License
[MIT License](LICENSE)

## Contact
For questions or collaboration, please reach out to the authors.


