<h1>
<table border="0">
  <tr>
    <td>
      <img src="images/SPECTRA_logo.png" width="80" height="80">
    </td>
    <td>
      <div align="center">
      <b>SPECTRA: Scalable Platform for Extracting</b><br>
      <b>Constraint-based Top-down Reconstructions and Analysis</b>
      </div>
    </td>
  </tr>
</table>
</h1>

## Overview

SPECTRA is a unified, scalable framework for reconstructing metabolic models from multi-omics data. It integrates diverse optimization strategies (linear and mixed-integer programming) to extract metabolic networks at multiple biological scales—from minimal reactomes and context-specific models (CSMs) to microbial community models and multi-tissue models. 


## Installation & Setup

### Prerequisites

- MATLAB
- [COBRA toolbox](https://opencobra.github.io/cobratoolbox/stable/index.html)
- CPLEX or Gurobi (for optimization)

### Installation

```bash
git clone https://github.com/NiravBhattLab/SPECTRA.git
cd SPECTRA
```

## Quick Start Guide
#### Step 1: Prepare The Inputs

The omics data can be used to derive the following three inputs to SPECTRA.
| Parameter | Purpose |
|---|---|
| **Core reactions** | These reactions will be present in the final model |
| **Weights** | Numbers that determine the relative importance of the reactions to be included (or excluded) |
| **Box contraints** | Lower and upper bounds to the reactions in the universal model|

The following code demonstrates the same for a toy network model:

```MATLAB
  % loading the universal model
  model = three_pathway_toy_model(); % A toy model with eight metabolites and eleven reactions

  % defining the box constraints 
  model.ub(:) = 10; 
  model.ub([1,6,9]) = 1; % constraining the media reactions

  % defining the core reactions
  core = 5; % the reaction indices of the core reactions

  % defining the weights of the reactions (For this example, we are using tradeOff as the network inference method)
  weights = [1,1,-2,-1,1,1,2,-1,-2,0,1]'; % This has to be in same size as the number of reactions in the universal model
```

#### Step 2: Choose The Network Inference Formulation

| Formulation | Objective | Use Case | Input Weights |
|---|---|---|---|
| **minNetLP** | Minimize total flux | Extracting minimal networks preserving core reactions | Non-negative weights |
| **minNetMILP** | Minimize reaction count |  Extracting minimal networks preserving core reactions | Non-negative weights |
| **minNetDC** | Cardinality minimization |  Extracting minimal networks preserving core reactions | Non-negative weights |
| **tradeOff** | Maximize weighted reactions | Balancing positive vs. negative evidence | Any real values |
| **optimBiomass** | Maximize biomass + minimize flux | Biological growth optimization | Non-negative weights |

#### Step 3: Model reconstruction

```matlab
% For flux consistent universal models
consType = 'stoichiometry'; % to use stoichiometric constraints or topology constraints
probType = 'tradeOff'; % network inference method
altSolMethod = {}; % method to obtain alternate solutions
nSol = 1; % number of alternate solutions required
solveTime = 7200; % maximum time for solving MILP problem
tol = 1e-5; % defining the small positive number
[Model_tradeOff] = spectraME(model,core,tol,consType,weights,nSol,altSolMethod,probType,solveTime); % model extraction

# For flux inconsistent universal models (simultaneous consistency checking + extraction)
consType = 'stoichiometry'; % to use stoichiometric constraints or topology constraints
probType = 'tradeOff'; % network inference method
altSolMethod = {}; % method to obtain alternate solutions
nSol = 1; % number of alternate solutions required
solveTime = 7200; % maximum time for solving MILP problem
tol = 1e-5; % defining the small positive number
[Model_tradeOff] = spectraCCME(model,core,tol,consType,weights,nSol,altSolMethod,probType,solveTime); % model extraction

```

#### Step 4: To Generate Alternative Solutions

```matlab
  % Pathway exclusion (for MILP formulations like tradeOff and minNetMILP)
  consType = 'stoichiometry'; % to use stoichiometric constraints or topology constraints
  probType = 'tradeOff'; % network inference method
  altSolMethod = 'pathwayExclusion'; % method to obtain alternate solutions
  nSol = 5; % number of alternate solutions required
  solveTime = 7200; % maximum time for solving MILP problem
  tol = 1e-5; % defining the small positive number
  [Model_tradeoff_altSols] = spectraME(model,core,tol,consType,weights,nSol,altSolMethod,probType,solveTime); % models extraction

```
