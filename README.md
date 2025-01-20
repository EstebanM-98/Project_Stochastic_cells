# Project: [On stochasticity in cell dynamics]

This repository contains a python jupyter notebook and a report of the solution to the stochastics cells problem from the Sethna Entropy, Order Parameters, and Complexity book. 

## Contents

- stocastic_cell.ipynb: In this notebook you will find the problem statement and the developed functions that lead to the solution. The following plots were found for the solution.

**Comparison between theoretical steady state and numerical solution to the equation describing the reactions of monomers.**
![defect](https://github.com/EstebanM-98/Project_Stochastic_cells/blob/a3a441dcae8be7c7966b14f6be8e9977f52f4bb7/Images/Nvst_continuum.png)

**Comparison of the theoretical solution when considering a finite number of monomers and the solution by Monte Carlo**
![defect1](https://github.com/EstebanM-98/Project_Stochastic_cells/blob/75daf145041a562f334ed62b03277039165f6bb8/Images/Nvst_stoc_comp.png)

**Mean value of the percentage error as a function of the number of monomers using the montecarlo algorithm for large time.**
![defect2](https://github.com/EstebanM-98/Project_Stochastic_cells/blob/75daf145041a562f334ed62b03277039165f6bb8/Images/Mean_percentage_error_vsN_without_many_realizatios.png)

**Comparison between theoretical solution and numerical solution using the montecarlo algorthm considering many realization and taking the average.**
![defect3](https://github.com/EstebanM-98/Project_Stochastic_cells/blob/75daf145041a562f334ed62b03277039165f6bb8/Images/Mvst_many_realizations_comparison.png)

**Mean value of the percentage error as a function of the number of monomers using the montecarlo algorithm for large time considering many realization and taking the average.**
![defect4](https://github.com/EstebanM-98/Project_Stochastic_cells/blob/75daf145041a562f334ed62b03277039165f6bb8/Images/Mean_percentage_error_vsN_with_many_realizatios.png)


## Prerequisites

To run this notebook, you need:

- Python 3.7.2 or higher.
- The following Python libraries:
  - `numpy`
  - `matplotlib`
  - `pandas`
  - `scipy`

## How to Run the Notebook

### Option 1: Use Google Colab (Recommended)
You can run the notebook in Google Colab without installing anything locally:

1. Click the button below:
   [![Open in Google Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1KiGS7GFFGm1o5oKJjWtOs8WhnTbzIhDy)
   
2. Once the notebook is open in Colab, run the cells directly.

### Option 2: Clone the Repository and Run Locally

1. Clone this repository:
   ```bash
   git clone https://github.com/your_username/your_repositor
