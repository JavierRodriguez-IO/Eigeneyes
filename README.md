# Eigeneyes

Code accompanying the paper *Eigeneyes: A statistical model of whole-eye geometry and its applications* (Vadillo, de la Peña, Rodríguez-Sánchez, de Castro, Marcos and Martínez-Enríquez).

MATLAB code to build a statistical shape model of the whole eye ("eigeneyes") from ocular biometry surfaces using principal component analysis (SVD), and to evaluate it with subject-grouped 10 × 10 cross-validation.

## Requirements

- MATLAB (tested with R2025b)
- Statistics and Machine Learning Toolbox (`cvpartition`)

## Repository structure

```
Code/
├── Eigeneyes_script_limpio.m          # main script
└── functions/
    ├── eigeneyes_image_registration.m # sets the coordinate origin
    └── eigeneyes_SVD_cartesian.m      # SVD / eigeneye computation
Data/                                  # input dataset (not included)
Results/                               # outputs are written here
```

## Inputs

`Data/eigeneyes_complete_dataset_all_right.mat`, containing:

| Variable | Description |
|---|---|
| `dataZ_eigeneyes` | Points × eyes matrix with the Z coordinates (mm) of the stacked ocular surfaces (anterior/posterior cornea, anterior/posterior lens, retina) |
| `dataX_eigeneyes`, `dataY_eigeneyes` | X and Y sampling coordinates of the grid |
| `Age` | Age of each eye (years) |
| `ID` | Subject identifier per eye (both eyes of a subject share the same ID) |

All eyes are expressed as right eyes (left eyes mirrored, x → −x).

The clinical data are not distributed with this repository because of patient privacy restrictions.

## Process

Run `Code/Eigeneyes_script_limpio.m`. It:

1. Loads the dataset and registers every eye with the origin at the retina (`eigeneyes_image_registration`).
2. Computes the final eigeneye model on all eyes: mean eye, principal components (eigeneyes), per-eye coefficients, cumulative explained variance and reconstruction RMSE with K = 10 eigeneyes (`eigeneyes_SVD_cartesian`).
3. Runs 10 × 10 cross-validation grouped by subject (both eyes of a subject always fall in the same fold). In each fold, the PCA is fitted on the training eyes only and the test eyes are reconstructed with K = 1–10 eigeneyes, recording the reconstruction RMSE and the explained variance (PVar).
4. Summarises the mean ± SD of RMSE and PVar for each K and plots them.

The random seed is fixed (`rng(42)`) for reproducibility.

## Outputs

Written to `Results/`:

| File | Content |
|---|---|
| `Eigeneyes_modelo_final.mat` | Final model: mean eye (`mn`), eigeneyes (`PC`), coefficients (`signals`), cumulative explained variance, reconstructed eyes and per-eye RMSE |
| `Eigeneyes_10x10_CV.mat` | Full cross-validation results (per repetition, fold, K and eye) |
| `Eigeneyes_10x10_CV_results.csv` | Summary table: RMSE (mm) and PVar (%) mean ± SD for K = 1–10 |
| `Eigeneyes_10x10_CV_RMSE_PVAR.png` | RMSE and PVar versus number of eigeneyes |

## License

Copyright (C) 2026 Estela Vadillo, Álvaro de la Peña, Javier Rodríguez-Sánchez, Alberto de Castro, Susana Marcos, Eduardo Martínez-Enríquez.

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. See [LICENSE](LICENSE) for the full text.
