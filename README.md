# Generalized propensity score and covariate matching for multilevel treatments in hierarchical data: Rebuilding Guo, Liu, and Wang (2022) in R

Qiu Wang, Syracuse University · v1.0 · November 2026

## Key findings

Guo, Liu and Wang (2022) proposed covariate matching (CM) and generalized propensity score (GPS) matching to estimate pairwise effects across multilevel treatment combinations, and reported simulations showing that both reduce the bias of the raw group difference.

- **Published vs rebuilt.** Their two simulations (4 and 9 treatment combinations) report that matching removes most of the raw-difference bias. A rebuild written from the paper's text alone, without the original code, confirms this: mean absolute bias falls from 0.6–0.9 (raw difference) to 0.2–0.3 (CM) and 0.03–0.10 (GPS).
- **Where the numbers differ.** The exact values do not match. Our GPS bias is much smaller than the paper's Table 1 (0.026 against 0.29 in Simulation 1), and in Simulation 2 our bias pattern tracks the published one with the opposite sign (correlation about −0.95 for the raw difference).
- **Gap partly explained.** The paper leaves thirteen points open (D1–D13). Reading Eq. (4) literally, as printed, closes about half of the GPS gap (mean absolute bias 0.026 to 0.159, against 0.293 published); the rest is unexplained, and we cannot say which unreported choice accounts for it.
- **The method stands.** The paper's estimators are implemented as described and are unchanged. Both remove most of the raw-difference bias; CM is more efficient (about half GPS's MSE) and GPS is less biased, a bias–variance trade-off.

## Why this replication

The paper offers estimators for treatment effects when students are nested in schools, but it provides no code or data and leaves several simulation details unspecified. This rebuild checks which reported results hold up, documents the choices needed to run it (D1–D13), and serves as a baseline for extending matching to clustered multilevel data.

## Background

Everything here was written from the published description, without the original code. The paper itself is not stored in this repository:

> Guo, S., Liu, J., & Wang, Q. (2022). Effective learning during COVID-19: Multilevel covariates matching and propensity score matching. *Annals of Data Science, 9*(5), 967–982. https://doi.org/10.1007/s40745-022-00392-x

Individual-level covariates X and school-level covariates Z are combined with treatments at two levels; every combination of the two treatment levels is one treatment group, so the problem becomes a comparison of A groups and all A(A − 1)/2 pairwise effects.

**How the paper handles the hierarchy.** The paper's setting is hierarchical, with students nested in schools and a treatment at each level, but its estimation is single-level:

- **Treatment.** The student-level and school-level treatments are combined into one multi-valued treatment, T = 1, …, A; Section 4 treats the problem as a multiple-treatment case.
- **Matching.** Both methods match across all units on student covariates X and school covariates Z. Cluster membership is not used: there are no cluster IDs, random effects or within-school matching. Without Z, the methods reduce to those of Yang et al. (2016) for multi-level treatments in non-hierarchical data.
- **Simulations.** Z is drawn separately for each unit, so the simulated data are not nested in schools.

This rebuild follows the paper's design. Matching that uses the clustering, such as matching within schools first, is left for future work.

## Design (as in Section 5)

- **Covariates.** X = (X1, …, X4) ~ multivariate normal with the paper's covariance matrix; Z1 ~ Uniform(−3, 3) (D2), Z2 ~ χ²(1), Z3 ~ Bernoulli(0.5).
- **Treatment.** Multinomial logit with η₁ = 0 and η_t = 0.1 × a sum of linear and squared covariate terms (D1). Simulation 1 has A = 4 groups of 100 units (6 effects); Simulation 2 has A = 9 groups of 200 units (36 effects).
- **Outcomes.** Y(t) is linear in (X, Z) with group-specific coefficients, plus N(0, 1) noise.
- **Estimators.** DIF, the raw difference in group means; CM, one Mahalanobis nearest neighbour with replacement on all 7 covariates; GPS, one nearest neighbour on the score from a linear multinomial logit (misspecified on purpose).
- **Inference.** Percentile bootstrap, B = 200, resampling within treatment groups and refitting the GPS model each time (D8).
- **Scale.** 500 replications per simulation. Bias, MSE and coverage are reported against two truths: the population ATE (τ\_pop) and the balanced-design ATE (τ\_bal) (D4, D9).

The full set-up is in [`docs/spec.md`](docs/spec.md).

## What differed

The paper leaves thirteen points open or inconsistent. The rebuild documents a choice for each in [`docs/spec.md`](docs/spec.md) (D1–D13). The ones that matter most:

1. **Treatment probabilities (D1).** The printed denominator, 1 + Σ exp(η), makes the probabilities sum to less than 1. We use the standard multinomial logit.
2. **Z1 (D2).** "N[−3, 3]" is read as Uniform(−3, 3).
3. **Fixed group sizes (D3, D9).** Units are drawn from a large pool by the multinomial probabilities and the first 100 or 200 per group are kept. This sample is not drawn from the population, so CM and GPS estimate the balanced-design ATE; both truths are reported.
4. **Covariance matrix for matching (D6).** The paper builds it from X only but matches on X and Z; we use all 7 covariates.
5. **Eq. (4) (D7, D13).** The printed formula uses the score for t in both terms. The main GPS uses the score for t′ in the second term; the literal reading is reported as a sensitivity check in [`docs/sensitivity.md`](docs/sensitivity.md).

## Results

Mean absolute bias (95% CI coverage), averaged over pairs. Per-pair values are in [results/summary.md](results/summary.md); the full comparison with the paper is in [docs/replication_verdict.md](docs/replication_verdict.md), and pair-by-pair tables against the paper's Tables 1 and 3 are in [docs/two_way.md](docs/two_way.md).

| Estimator | Sim 1, vs τ\_pop | Sim 1, vs τ\_bal | Sim 2, vs τ\_pop | Sim 2, vs τ\_bal |
| --- | --- | --- | --- | --- |
| DIF | 0.588 (85%) | 0.612 (84%) | 0.941 (78%) | 0.927 (79%) |
| CM | 0.187 (91%) | 0.223 (90%) | 0.295 (78%) | 0.262 (81%) |
| GPS | 0.026 (100%) | 0.072 (100%) | 0.099 (100%) | 0.096 (100%) |

- **Efficiency.** GPS has about twice CM's MSE (against τ\_pop: 0.68 vs 0.36 in Simulation 1, 0.63 vs 0.29 in Simulation 2), and its intervals are about 30% too wide, which gives its 100% coverage.
- **Against the paper.** Simulation 1 GPS bias averages 0.026 here and 0.29 in the paper's Table 1. In Simulation 2 our bias tracks Table 3 with the opposite sign (correlations −0.95, −0.76 and −0.87 for DIF, CM and GPS).
- **Sensitivity (D13).** Reading Eq. (4) literally raises the Simulation 1 GPS bias to 0.159, about half the way to the paper's value.

## Reproduce

From the repository root (R 4.5; packages MASS and nnet):

| Step | Command | Time |
|---|---|---|
| 0. True effects (needed by the summaries) | `Rscript R/check_truth.R` | about 1 min |
| (optional) Tests | `Rscript R/test_estimators.R` (`quick` for unit tests only) | about 6 min |
| 1. Simulations 1 and 2, then summaries | `Rscript R/run_all.R` (options: `R=500 B=200 workers=10`) | On 10 workers: Simulation 1 about 15 min, Simulation 2 about 4 h 48 min, full run about 5 h (log: 14m28s, 4h48m19s, 5h03m08s) |
| 2. Rebuild summaries only | `Rscript R/summarize.R` | seconds |
| 3. Figures and Simulation 2 comparison | `Rscript R/plots.R` | seconds |
| 4. D13 sensitivity (Simulation 1) | `Rscript R/run_gps_literal.R`, then `Rscript R/run_gps_literal2.R` (`check` runs the self-tests) | |

Each replication is saved to `results/reps/simK/rep_###.rds` when it finishes. Re-running `run_all.R` skips replications that are already done, so an interrupted run resumes where it stopped. It refuses to mix checkpoints made with different settings (R, B, seed, optimizer tolerance). Progress goes to standard output; redirect it to `results/run_log.txt` for a background run.

Each replication uses its own pre-generated L'Ecuyer-CMRG random-number stream, so results do not depend on the number of workers or the order in which replications run. The scripts use only base R, `parallel` and `nnet`, plus `MASS` for data generation. The R version and session details of the full run are in `results/run_session_info.txt`; the `nnet` and `MASS` versions used (7.3.20 and 7.3.65) are on the second line of `results/run_log.txt`. No `renv` lockfile is included.

## Files

```
README.md, LICENSE, CITATION.cff
docs/spec.md                  simulation specification and decisions D1-D13
docs/replication_verdict.md   results and decisions D1-D13
docs/two_way.md               our results against the paper's published tables
docs/sensitivity.md           D13 sensitivity check (literal Eq. (4))
R/                            scripts (data generation, estimators, runners, summaries)
results/                      output tables, figures and summary.md (large .rds files are git-ignored)
```

| Path | Contents |
|---|---|
| [docs/spec.md](docs/spec.md) | simulation specification, with the decisions D1–D13 |
| [docs/replication_verdict.md](docs/replication_verdict.md) | results and decisions D1–D13: what replicates, what does not, and what the rebuild cannot show |
| [docs/two_way.md](docs/two_way.md) | two-way comparison: our results and the paper's published Tables 1 and 3, pair by pair, with the overlay figures |
| [docs/sensitivity.md](docs/sensitivity.md) | D13 sensitivity check: GPS_lit and GPS_lit2 (Eq. (4) read literally), Simulation 1 only |
| [results/summary.md](results/summary.md) | compact per-pair tables reformatted from the summary CSVs: both simulations against τ\_pop and τ\_bal, and against the paper's Tables 1 and 3 |
| [results/true_effects_sim1.csv](results/true_effects_sim1.csv), [results/true_effects_sim2.csv](results/true_effects_sim2.csv) | tau_pop and tau_bal for every pair, with Monte Carlo SEs |
| [results/sim1_summary.csv](results/sim1_summary.csv), [results/sim2_summary.csv](results/sim2_summary.csv) | bias, absolute bias, SD, MSE, coverage and CI width for each estimator, pair and target (`pop` or `bal`) |
| [results/sim1_diagnostics.csv](results/sim1_diagnostics.csv), [results/sim2_diagnostics.csv](results/sim2_diagnostics.csv) | per replication: run time, GPS convergence, matching distances, share of p_hat < 0.01 |
| [results/sim1_vs_paper.csv](results/sim1_vs_paper.csv), [results/sim2_vs_paper.csv](results/sim2_vs_paper.csv) | our numbers next to the paper's Table 1 and Table 3 |
| [results/fig1_bias_sim2_pop.png](results/fig1_bias_sim2_pop.png), [results/fig1_bias_sim2_bal.png](results/fig1_bias_sim2_bal.png) | Fig. 1 replica: signed bias of the 36 Simulation 2 effects |
| [results/fig1_bias_sim2_pop_vs_paper.png](results/fig1_bias_sim2_pop_vs_paper.png), [results/fig1_bias_sim2_pop_vs_paper_flipped.png](results/fig1_bias_sim2_pop_vs_paper_flipped.png) | the same, with the paper's values as dashed lines: as published, and multiplied by -1 (see "Comparison with the paper" in [docs/spec.md](docs/spec.md)) |
| [results/sim2_bias_correlations.csv](results/sim2_bias_correlations.csv) | correlation of our bias with the paper's across the 36 pairs, paper sign as published and flipped |
| [results/fig2_table5_redrawn.png](results/fig2_table5_redrawn.png) | the paper's Fig. 2 redrawn from its Table 5 (not replicated) |
| [results/paper_table1.csv](results/paper_table1.csv), [results/paper_table3.csv](results/paper_table3.csv), [results/paper_table5.csv](results/paper_table5.csv) | values transcribed from the paper's published tables, for comparison only |
| [results/run_log.txt](results/run_log.txt), [results/run_session_info.txt](results/run_session_info.txt) | log and `sessionInfo()` of the full run |
| results/reps/ (local only, not in the repo) | per-replication checkpoints |
| results/reps_gps_literal/, results/reps_gps_literal2/ (local only, not in the repo) | per-replication checkpoints of the D13 GPS_lit and GPS_lit2 runs |
| [results/gps_literal_sim1.csv](results/gps_literal_sim1.csv) | D13 sensitivity, Simulation 1: GPS vs GPS_lit (Eq. (4) read literally for units outside groups t and t') per pair, with the paper's Table 1 GPS values |
| [results/gps_literal2_sim1.csv](results/gps_literal2_sim1.csv) | D13 sensitivity, Simulation 1: GPS vs GPS_lit vs GPS_lit2 (fully literal Eq. (4), also for units in group t) per pair, with the paper's Table 1 GPS values |
| [R/gps_literal.R](R/gps_literal.R) | the GPS_lit estimator, reproducing the main run's random-number sequence |
| [R/gps_literal2.R](R/gps_literal2.R) | the GPS_lit2 estimator (requires [R/gps_literal.R](R/gps_literal.R)) |
| [R/run_gps_literal.R](R/run_gps_literal.R) | self-tests (`check`) and the paired GPS_lit run on the main Simulation 1 datasets; writes [gps_literal_sim1.csv](results/gps_literal_sim1.csv) |
| [R/run_gps_literal2.R](R/run_gps_literal2.R) | self-tests (`check`) and the paired GPS_lit2 run; writes [gps_literal2_sim1.csv](results/gps_literal2_sim1.csv) |

The paper's PDF is not included. The `paper_table*.csv` files hold values transcribed from its published tables, for comparison only.

## Lesson

The GPS model's fitting tolerance matters when matching on fitted probabilities. With `nnet::multinom`'s defaults (reltol = 1e-8, maxit = 100), the fitted probabilities were not precise enough for the very small matching distances, and the GPS estimate depended on the optimizer's starting point. Fitting at reltol = 1e-14 with maxit = 10000 (D5) removed this, and a test now checks it.

## Conclusions

In these simulations, covariate matching and GPS matching both remove most of the bias in the raw group difference, as the paper reports. GPS had the smallest bias and CM the smallest mean squared error, so the better choice depends on whether bias or variance matters more. The published numbers cannot be reproduced exactly from the paper alone, because thirteen points are left open; reading Eq. (4) literally accounts for about half of the GPS discrepancy. These results are a baseline for extending the comparison to clustered, multilevel data.

## Acknowledgment

This is a re-implementation, written from the published description of Guo, Liu & Wang (2022). The original authors' code was not used.

Disclosure: the author of this repository is a coauthor of Guo, Liu and Wang (2022).

## AI assistance

The code was written with assistance from Claude Code (Anthropic).

## How to cite

Wang, Q. (2026). *Generalized propensity score and covariate matching for multilevel treatments in hierarchical data: Rebuilding Guo, Liu, and Wang (2022) in R* (Version 1.0) [Computer software]. https://github.com/wangqiu-syr/glw2022-matching-replication

DOI: to be added after the Zenodo release. Citation metadata are in [`CITATION.cff`](CITATION.cff); please also cite the original paper.

## License

MIT (see `LICENSE`).
