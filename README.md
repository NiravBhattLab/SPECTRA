<h1 align="center">
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

SPECTRA is a unified, scalable framework for reconstructing metabolic models from multi-omics data. It integrates diverse optimization strategies (linear and mixed-integer programming) to extract metabolic networks at multiple biological scales-from minimal reactomes and context-specific models (CSMs) to microbial community models and multi-tissue models. 

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

## Reproducing Manuscript Results

This repository contains the code necessary to reproduce the results and figures presented in the main manuscript and supplementary information. The scripts are organized into corresponding folders.

### 1. Build Cancer Models (Context-Specific Models)
Located in `BuildCancerModels/`
Build context-specific models using DepMap transcriptomics data.
- **Workflow**:
  1. `DepMapModels_Recon/PrepareDepMapData.m`: Extract transcriptomics data.
  2. `DepMapModels_Recon/modGenesIDConv.m`: Convert gene IDs in `UpdatedRecon3D.mat`.
  3. `DepMapModels_Recon/GiniReactionImportance.m`: Obtain core reactions from transcriptomics.
  4. `DepMapModels_Recon/BuildDepMapModels.m`: Build models using fastcore, swiftcore, and sprintcore methods.
- **Analysis**:
  - `coverage_core/`: Scripts like `Cancerrxnscoverage.m` to calculate cancer generic core reaction coverage.
  - `flux_results/` and `fva/`: Evaluate FVA span ratios and flux enrichment.

### 2. Runtime Comparisons
Located in `RuntimeComparisons/`
Use these scripts to benchmark the runtime performance of SPECTRA formulations.
- `checkModelQuality.m`: Helper script to evaluate outputs.
- `ForSPECTRA_CCResults.m` & `ForSPECTRA_MEResults.m`: Recreate speed metrics for Consistency Checker and Model Extraction respectively.

### 3. Gap-filling & Random Reaction Removal (AGORA)
Located in `RandomReactionRemovalAGORA/` and `gapFilledAGORA_DB/`
Evaluates the gap-filling performance and stability of SPECTRA using microbial models.
- `randomReactionRemovalResults.m` (in `Visualisations/`): Generates F1-score/Accuracy plots by continuously dropping required reactions.
- `GapFillAGORADraftModels.m`: Gapfilling experiments across varied model degradation setups.
- `gapFilledAGORA_DB/BuildUmodelFromAGORAdb.m`: Generates universal models for draft AGORA models. 

### 4. Deep Learning EC Match Verification
Located in `Verify_new_reactions_AGORA_using_DL/`
Verify the biochemical viability of gap-filled reactions.
- `getECmatch.m`: Matches deep-learning-based Enzyme Commission (EC) numbers with newly predicted reactions.

### 5. Linear Programming Consistency Checks (LP7 vs LPf)
Located in `LP7_vs_LPf/`
Evaluates fast implementations for locating unblocked/consistent reactions.
- `fc_lp7.m` & `LP7_LPf.m`: Benchmarking and demonstrating the robust efficiency of SPECTRA LP-based fast consistency checker vs existing tools.

### 6. Minimal Microbiome Generation
Located in `MinimalMicrobiome/` and `hCOM_community/`
- Scripts like `createCommModel.m` and `getMinMicrobiome.m` assemble multi-species microbial environments and prune them using optimization protocols to create minimal functionally-stable communities.

### 7. Visualisations & Figures
Located in `Visualisations/`
Contains Python Notebooks (`.ipynb`) and MATLAB scripts that were used to plot the final paper figures. 
- *Note:* Use `conda activate momi` before running `.ipynb` files like `Runtime_SPECTRA_CC.ipynb`, `tSNE_plot_rxns.ipynb`, etc.

---

## Quick Start Guide

### Step 1: Prepare The Inputs

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
  weights = [1,1,-2,-1,1,1,2,-1,-2,0,1]'; % Same size as the number of reactions
```

### Step 2: Choose The Network Inference Formulation

| Formulation | Objective | Use Case | Input Weights |
|---|---|---|---|
| **minNetLP** | Minimize total flux | Extracting minimal networks preserving core reactions | Non-negative weights |
| **minNetMILP** | Minimize reaction count |  Extracting minimal networks preserving core reactions | Non-negative weights |
| **minNetDC** | Cardinality minimization |  Extracting minimal networks preserving core reactions | Non-negative weights |
| **tradeOff** | Maximize weighted reactions | Balancing positive vs. negative evidence | Any real values |
| **optimBiomass** | Maximize biomass + minimize flux | Biological growth optimization | Non-negative weights |

### Step 3: Model reconstruction

```matlab
% For flux consistent universal models
consType = 'stoichiometry'; % to use stoichiometric constraints or topology constraints
probType = 'tradeOff'; % network inference method
altSolMethod = {}; % method to obtain alternate solutions
nSol = 1; % number of alternate solutions required
solveTime = 7200; % maximum time for solving MILP problem
tol = 1e-5; % defining the small positive number
[Model_tradeOff] = spectraME(model,core,tol,consType,weights,nSol,altSolMethod,probType,solveTime); 

% For flux inconsistent universal models (simultaneous consistency checking + extraction)
[Model_tradeOff] = spectraCCME(model,core,tol,consType,weights,nSol,altSolMethod,probType,solveTime); 
```

### Step 4: To Generate Alternative Solutions

```matlab
  % Pathway exclusion (for MILP formulations like tradeOff and minNetMILP)
  altSolMethod = 'pathwayExclusion'; % method to obtain alternate solutions
  nSol = 5; % number of alternate solutions required
  [Model_tradeoff_altSols] = spectraME(model,core,tol,consType,weights,nSol,altSolMethod,probType,solveTime); 
```
