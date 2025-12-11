---
---
---

# Progress Report 02 — Full Monte Carlo Engine + Preliminary Findings

**Project:** Monte Carlo evaluation of sample-selection estimators under distributional misspecification\
**Author:** Zeynep Bozoklu\
**Status date:** 2025-12-11

This note documents (i) completion of the Monte Carlo engine over the full scenario grid and (ii) first summary statistics on estimator performance, with a focus on sign reliability of the main slope coefficient.

------------------------------------------------------------------------

## 1. Implementation status

**Design & calibration (as in Report 01)**\
- 48 scenarios defined by\
- $n \in \{500, 2000\}$\
- target selection rate $p \in \{0.3, 0.6\}$\
- correlation $\rho \in \{0, 0.3, 0.6\}$ between outcome and selection errors\
- error family $\in \{\text{Normal}, t_3, t_5, \text{Mixture}\}$.\
- For each scenario, the selection intercept $\gamma_0$ is calibrated so that $\Pr(S = 1) \approx p$ via root finding on\
$\mathbb{E}[F_v(\gamma_0 + Z\gamma_{\text{slopes}})] - p = 0$.

**DGP (one dataset)**\
For a given scenario: - Regressors: $X = (x_1, x_2) \sim N(0, \Sigma_X)$ with $\text{corr}(x_1,x_2) = 0.3$.\
- Selection covariates: $Z = (1, x_1, x_2)$.\
- Errors $(u,v)$ drawn with correlation $\rho$ and family specified by the scenario\
(Normal, $t_3$, $t_5$; Mixture modifies the outcome error $u$).\
- Latent outcome and selection: $$
  Y^* = X\beta + u, \quad \beta = (1,\,-0.5);
  \qquad S = 1\{Z\gamma + v > 0\}.
  $$ - Observed outcome: $Y = Y^*$ if $S=1$, else $Y = \text{NA}$.

**Estimators implemented** 1. **Selected-sample OLS**\
- OLS of $Y$ on $(1,X)$ using only observations with $S=1$. 2. **Zero-imputation OLS**\
- Replace missing outcomes with 0 for $S=0$, then OLS of $Y$ on $(1,X)$ using the full sample. 3. **Heckman two-step**\
- Probit: $S$ on $Z$ (no extra intercept); compute $\hat z_i = Z_i\hat\gamma$.\
- Inverse Mills ratio: $\hat\lambda_i = \phi(\hat z_i)/\Phi(\hat z_i)$.\
- Outcome regression on selected sample: $Y$ on $(1,X,\hat\lambda)$.\
- Keep only the coefficients on $(1,X)$ for comparison.

**Monte Carlo engine** - For a single scenario: - Repeat R times: simulate one dataset, run each estimator, and store the coefficient vector $(\alpha,\beta_1,\beta_2)$. - True parameter vector used for bias/RMSE: $(0,1,-0.5)$. - Summary for each estimator: - mean, bias, RMSE per coefficient, - sign-failure rate for $\beta_1$ (fraction of replications where $\hat\beta_1$ has the wrong sign). - Full MC: - For this first run, use $R = 50$ replications per scenario. - Loop over all 48 scenarios and all three estimators. - Export tidy summary table to `output/results/mc_summary.rds` and `.csv` with one row per (scenario, estimator, coefficient).

------------------------------------------------------------------------

## 2. Preliminary findings (R = 50 per scenario)

Let $\beta_1$ be the coefficient on $x_1$ (the main “treatment/effect” of interest).\
Average sign-failure rates for $\beta_1$, aggregated across all $n, p, \rho$ within each error family:

-   **Selected-sample OLS**
    -   Normal: 0\
    -   $t_3$: 0\
    -   $t_5$: 0\
    -   Mixture: 0
-   **Zero-imputation OLS**
    -   Normal: 0\
    -   $t_3$: 0\
    -   $t_5$: 0\
    -   Mixture: 0
-   **Heckman two-step**
    -   Normal: ≈ 1.8% sign failures\
    -   $t_5$: ≈ 4.7%\
    -   Mixture: ≈ 3.5%\
    -   $t_3$: ≈ **9.8%**

So in this pilot design, both OLS variants are biased but **directionally very stable** (never flipping the sign of $\beta_1$), while Heckman is more accurate on average in benign settings but shows **nontrivial sign-instability**, especially under heavy tails.

------------------------------------------------------------------------

## 3. Illustration: “harsh” $t_3$ scenario

One particularly revealing scenario is:

-   $n = 500$, $p_{\text{select}} = 0.3$, $\rho = 0.6$, error family $t_3$.

For this scenario (ID S021), the Monte Carlo summaries for $\beta_1$ are:

-   **Selected-sample OLS**:
    -   bias ≈ −0.39 (downward bias from 1), RMSE ≈ 0.45, sign_fail = 0.
-   **Zero-imputation OLS**:
    -   bias ≈ −0.57, RMSE ≈ 0.58, sign_fail = 0.
-   **Heckman two-step**:
    -   bias ≈ +0.24, RMSE ≈ 1.81, **sign_fail ≈ 0.20** (≈20% of replications wrong sign).

This scenario illustrates the emerging trade-off: - Heckman can be closer to the truth on average but **high-variance and fragile**, occasionally delivering extreme and even sign-reversing estimates in non-normal, small-sample, strong-selection settings. - Simpler OLS-based corrections remain systematically biased but **never reverse the sign** of the effect in this design.

------------------------------------------------------------------------

## 4. Next steps

1.  Possibly refine the mixture specification depending on feedback from Vedant (whether mixture affects only $u$ or both $u$ and $v$).\
2.  Increase $R$ t
