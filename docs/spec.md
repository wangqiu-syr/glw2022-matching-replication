# Simulation spec: multilevel covariate matching vs. GPS matching

Source: Guo, Liu & Wang (2022), *Annals of Data Science* 9(5):967-982, https://doi.org/10.1007/s40745-022-00392-x. This file restates the set-up in our own words and fixes the points the paper leaves open (marked **D1-D8**, listed at the end). Defaults can be changed; confirm them before the full run.

## Goal

Replicate the two simulation studies (paper Section 5; Tables 1 and 3) in R and compare three estimators of the pairwise average treatment effect across treatment combinations:

- **DIF**: raw difference in group means
- **CM**: covariate matching
- **GPS**: generalized propensity score matching

## Notation

- X = (X1, X2, X3, X4): individual-level covariates. Z = (Z1, Z2, Z3): school-level covariates.
- T in {1, ..., A}: treatment combination. A = 4 in Simulation 1, A = 9 in Simulation 2.
- Estimand: tau(t, t') = E[Y(t) - Y(t')] for every pair t < t'. This gives 6 effects in Simulation 1 and 36 in Simulation 2.

## Covariates (both simulations)

X ~ MVN(0, Sigma), with

```
Sigma = [  2.0   1.0  -1.0   0.5
           1.0   1.0  -0.5   0.2
          -1.0  -0.5   1.0  -0.3
           0.5   0.2  -0.3   2.0 ]
```

(positive definite; smallest eigenvalue about 0.31). Independent of X and of each other:

- Z1 ~ Uniform(-3, 3)  **(D2)**
- Z2 ~ chi-square(1)
- Z3 ~ Bernoulli(0.5)

## Treatment assignment

Multinomial logit, with eta_1 = 0 as the baseline:

p_t(x, z) = exp(eta_t) / sum over s of exp(eta_s)  **(D1)**

**Simulation 1 (A = 4):**

| t | eta_t / 0.1 |
|---|---|
| 1 | 0 (eta_1 = 0) |
| 2 | x1^2 + x2 + x3 + x4 + z1 + z2 + z3 |
| 3 | x1 + x2^2 + x3^2 + x4 + z1 + z2 + z3 |
| 4 | x1 + x2^2 + x3^2 + x4^2 + z1 + z2 + z3 |

**Simulation 2 (A = 9):** eta_1 = 0, and eta_t = 0.1 times:

| t | expression |
|---|---|
| 2 | x1^2 + x2 + x3 + x4 + z1 + z2 + z3 |
| 3 | x1 + x2^2 + x3^2 + x4 + z1 + z2 + z3 |
| 4 | x1 + x2^2 + x3^2 + x4^2 + z1 + z2 + z3 |
| 5 | x1 + x2^2 + x3^2 + x4^2 + z1^2 + z2 + z3 |
| 6 | x1 + x2^2 + x3^2 + x4^2 + z1^2 + z2^2 + z3 |
| 7 | x1 + x2^2 + x3^2 + x4^2 + z1^2 + z2^2 + z3^2 |
| 8 | x1^2 + x2^2 + x3^2 + x4^2 + z1^2 + z2^2 + z3^2 |
| 9 | x1 + x2 + x3 + x4 + z1 + z2 + z3 |

Only the index t matters. Ignore the cell labels (Aa, Ba, ...) in the paper's Table 2, which contain numbering typos.

## Potential outcomes

epsilon ~ N(0, 1), drawn independently for each unit. Observed outcome Y = Y(T).

**Simulation 1:** Y(t) = beta_t' (x1, x2, x3, x4, z1, z2, z3) + epsilon, with

| t | beta_t |
|---|---|
| 1 | (1, 1, 1, 1, 1, 1, 1) |
| 2 | (2, 3, 1, 2, 2, 2, 2) |
| 3 | (3, 1, 2, 1, 3, 1, 3) |
| 4 | (4, 1, 3, 4, 4, 1, 4) |

**Simulation 2:** Y(t) = c_t' (x1, x2, x3, x4, z1, z2 - 1, z3 - 0.5) + epsilon, with

| t | c_t |
|---|---|
| 1 | (1, 1, 1, 1, 1, 1, 1) |
| 2 | (2, 3, 1, 2, 2, 2, 2) |
| 3 | (3, 1, 2, 1, 3, 1, 3) |
| 4 | (4, 1, 3, 4, 4, 1, 4) |
| 5 | (5, 1, 1, 1, 1, 1, 1) |
| 6 | (6, 3, 1, 2, 2, 2, 9) |
| 7 | (7, 1, 8, 1, 3, 1, 8) |
| 8 | (8, 1, 5, 2.5, 4, 1, 7) |
| 9 | (9, 1, 7, 1.5, 4, 1, 7) |

The paper prints no epsilon in Y(1) for Simulation 2. Treat that as a typo and add it.

## Sample design  (D3)

| | Simulation 1 | Simulation 2 |
|---|---|---|
| Units per treatment combination | 100 | 200 |
| Total n | 400 | 1800 |
| Replications | 500 | 500 |

The paper fixes the group sizes but also assigns treatment by a multinomial model. Default: draw a large pool of units, assign T from p_t, then keep the first n_g units in each group (pool size large enough that every group fills; check this in code).

## True effects  (D4)

Target = population ATE over the pooled population, tau(t, t') = (beta_t - beta_t')' mu, where mu = E[(X, Z)] (use the centred version in Simulation 2). With Z1 mean 0:

- **Simulation 1:** mu = (0, 0, 0, 0, 0, 1, 0.5), so tau(1, 2) = -1.5, and so on.
- **Simulation 2:** all 36 true effects are 0, because z2 and z3 are centred in the outcome equations.

Verify these closed forms against a Monte Carlo average on a very large population (about 1e6 units) before trusting them.

## Estimators

For each replication and each pair (t, t'):

**DIF.** mean(Y | T = t) minus mean(Y | T = t').

**CM.** For every unit i in the sample and every level s, impute Y_i(s) = Y_i if T_i = s. Otherwise use the Y of the nearest unit j with T_j = s, where distance is Mahalanobis on all 7 covariates (X and Z together, V = pooled sample covariance)  **(D6)**. Use 1 match, with replacement, and break ties at random. Then tau_hat_CM(t, t') = mean over all i of [Y_i(t) - Y_i(t')].

**GPS.** Fit a multinomial logit of T on x1..x4, z1, z2, z3 (linear terms only, e.g. `nnet::multinom`)  **(D5)**. This is deliberately misspecified, because the true eta contains squares. For every unit i and level s, impute Y_i(s) = Y_i if T_i = s. Otherwise use the Y of the unit j with T_j = s that minimizes |p_hat_s(x_j, z_j) - p_hat_s(x_i, z_i)|. Then tau_hat_GPS(t, t') = mean over all i of [Y_i(t) - Y_i(t')]. Use p_hat_t' for the t' term. The paper's Eq. (4) prints p_t in both terms, which looks like a typo  **(D7)**.

## Inference and coverage  (D8)

The paper does not say how its confidence intervals were built. Default: nonparametric bootstrap, percentile 95% CI, B = 200 for the final run (B = 20 for tests). Note that bootstrapping matching estimators is not formally valid (Abadie and Imbens, 2008), so coverage may differ from the paper. Sensitivity option: Abadie-Imbens standard errors for CM.

## Outputs

Per estimator and pair, over the replications:

- signed bias = mean(tau_hat) - tau
- absolute bias = |mean(tau_hat) - tau| (Table 1 reports this)
- MSE
- coverage of the 95% CI (%)

Write tidy summaries to `results/sim1_summary.csv` and `results/sim2_summary.csv`. Write per-replication estimates to `results/*.rds` (git-ignored). Reproduce the bias plot of the 36 effects (paper Fig. 1) for Simulation 2.

## Implementation notes

- Folder plan: `R/dgp.R`, `R/estimators.R`, `R/run_sim1.R`, `R/run_sim2.R`, `R/summarize.R`.
- Per-replication seeds, parallelized with `future.apply` or `parallel`. Test with 20 replications before 500. Record package versions with `renv`.
- Sanity checks: p_t sums to 1; every group reaches its target size; closed-form tau matches Monte Carlo; GPS overlap looks sensible.
- Commit after each stage.

## Out of scope

Paper Section 6 (the learning-mode application). Its data-generating process is not described, so Table 5 cannot be reproduced. For reference, the cell sizes in the paper's Table 4 are 8, 56, 118, 93, 75 and 70.

## Open decisions and defaults

| ID | Issue | Default |
|---|---|---|
| D1 | Paper's denominator is 1 + sum(exp(eta)), which makes probabilities sum below 1 | Standard multinomial logit, no extra 1 |
| D2 | "Z1 ~ N[-3, 3]" is ambiguous | Uniform(-3, 3). Its mean of 0 is consistent with the centred z2, z3 terms, which make all true effects 0 in Simulation 2 |
| D3 | Fixed group size versus random assignment | Large pool, assign by p_t, keep n_g per group |
| D4 | True tau not stated | Pooled-population ATE, closed form checked by Monte Carlo |
| D5 | GPS model not stated | Linear multinomial logit (misspecified on purpose) |
| D6 | Paper defines V from X only but matches on X and Z | Mahalanobis on (X, Z) jointly |
| D7 | Eq. (4) uses p_t in both terms | Use p_t' for the t' term |
| D8 | CI method not stated | Percentile bootstrap, B = 200 |

Exact agreement with the published tables should not be expected, because of these open points. The qualitative pattern should hold: DIF is clearly biased, while CM and GPS reduce the bias substantially.
