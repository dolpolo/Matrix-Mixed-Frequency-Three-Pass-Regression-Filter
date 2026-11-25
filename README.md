# Matrix MF-TPRF

Replication code for **“Matrix Three-Pass Regression Filter: A Nowcasting
Framework for Mixed-Frequency Multidimensional Data.”**

## What it does

The Matrix Mixed-Frequency Three-Pass Regression Filter (Matrix MF-TPRF)
nowcasts quarterly targets from monthly predictors arranged as a matrix
(units × indicators). It uses quarterly proxies to extract target-relevant
factors while preserving the matrix structure, then combines them with
temporal aggregation and U-MIDAS regressions to handle mixed frequencies and
produce unit-specific real-time nowcasts.

## Run

Requires R 4.3 or later. From the repository root:

```bash
cd Matrix-MF-TPRF_Replication
Rscript install_dependencies.R
Rscript run_all.R
```

The package contains the source code, required inputs, and replication
results. See [its README](Matrix-MF-TPRF_Replication/README.md) for individual
workflows and further details.
