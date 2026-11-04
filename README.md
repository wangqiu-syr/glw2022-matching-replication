# ClaudeCodeTry

R simulation of two multilevel-treatment estimators for observational data:

- **CM**: covariate matching (Mahalanobis nearest neighbor)
- **GPS**: generalized propensity score matching

The set-up follows the simulation studies in Guo, Liu & Wang (2022), *Annals of Data Science*, 9(5), 967-982, https://doi.org/10.1007/s40745-022-00392-x. The paper itself is not stored in this repository.

## Status

Set-up stage. The full specification is in [`docs/spec.md`](docs/spec.md), including the open decisions that still need to be confirmed.

## Layout

```
README.md
docs/spec.md      simulation specification
R/                scripts (data generation, estimators, runners, summaries)
results/          output tables and figures (large .rds files are git-ignored)
```

## How to run

To be filled in once the scripts exist. Planned entry points: `R/run_sim1.R` and `R/run_sim2.R`.

## Reproducibility

Random seeds are set per replication. Package versions will be recorded with `renv`.
