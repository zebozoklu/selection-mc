# Progress Report 01 — Monte Carlo Design + Calibration (No Estimation Code Yet)

**Project:** Monte Carlo evaluation of sample-selection estimators under distributional misspecification\
**Author:** Zeynep Bozoklu\
**Status date:** 2025-12-11\
**Scope of this report:** Design objects + selection-rate calibration only. No DGP simulation loop or estimators implemented yet.

------------------------------------------------------------------------

## 1. Objective

Build a Monte Carlo framework to compare sample-selection estimators—especially the Heckman two-step—when the joint normality assumption fails (heavy tails / outliers / skew-type deviations). The focus is finite-sample performance and sign reliability of a target coefficient.

------------------------------------------------------------------------

## 2. Implementation blueprint (layered structure)

The codebase is organized conceptually into four layers:

1)  **Scenario / Design layer**\

-   **Input:** vectors of design knobs $(n,\; p_{\text{select}},\; \rho,\; \text{error family})$\
-   **Output:** scenario grid table (one row per scenario, unique `scenario_id`)

2)  **Calibration layer (implemented)**\

-   **Input:** scenario row + $\gamma_{\text{slopes}}$ + large draw of $Z$\
-   **Output:** calibrated selection intercept $\gamma_0$ such that $\Pr(S=1)\approx p_{\text{select}}$

3)  **DGP layer (next)**\

-   **Model:**\
    $$
    Y^* = X\beta + u,\quad
    S = 1\{Z\gamma + v > 0\},\quad
    Y = Y^* \text{ observed iff } S=1.
    $$
-   **Output:** one simulated dataset $(Y,S,X,Z)$

4)  **Estimator + Monte Carlo + Reporting layers (later)**\

-   Estimators: selected-sample OLS, zero-imputation OLS, Heckman 2-step (optional series correction)\
-   Monte Carlo engine: repeat simulate→estimate across replications/scenarios\
-   Reporting: bias/RMSE/sign-failure tables + plots

------------------------------------------------------------------------

## 3. Scenario grid (implemented)

A baseline design grid is constructed with:

-   Sample size: $n \in \{500, 2000\}$
-   Target selection rate: $p_{\text{select}} \in \{0.3, 0.6\}$
-   Dependence between unobservables: $\rho \in \{0, 0.3, 0.6\}$
    -   $\rho=0$ serves as a key sanity-check case (selection exists but not on unobservables).
-   Error family: $\{\text{normal}, t_3, t_5, \text{mixture}\}$

This yields **48 scenarios** total (balanced: 12 per error family).\
The scenario table stores needed distribution parameters (e.g., `df` for $t$, mixture parameters such as mixing probability/shift).

------------------------------------------------------------------------

## 4. Selection-rate calibration via $\gamma_0$ (implemented + tested)

### 4.1 Calibration target

For each scenario, choose $\gamma_0$ (the intercept in the selection index) so that: $$
\Pr(S=1) \approx p_{\text{select}}.
$$

### 4.2 Method

Using: $$
\Pr(S=1) = \mathbb{E}\big[F_v(\gamma_0 + Z\gamma_{\text{slopes}})\big],
$$ the expectation is approximated with a large simulated draw of $Z$ (e.g., 200,000 draws).\
Given the scenario’s selection-error CDF $F_v$ (normal or $t_\nu$), $\gamma_0$ is solved by 1D root-finding (`uniroot`), leveraging monotonicity in $\gamma_0$.

### 4.3 Sanity-check results (illustrative)

With $Z=(x_1,x_2)\sim N(0,1)$ and $\gamma_{\text{slopes}}=(0.8,-0.4)$, calibration produced:

-   **Normal** selection error $v$:
    -   $p=0.3 \Rightarrow \gamma_0 \approx -0.704$\
    -   $p=0.6 \Rightarrow \gamma_0 \approx 0.340$
-   **Heavy tails** in $v$:
    -   $t_3:\ \gamma_0 \approx -0.792$ (p=0.3), $0.380$ (p=0.6)\
    -   $t_5:\ \gamma_0 \approx -0.755$ (p=0.3), $0.363$ (p=0.6)

**Interpretation:** $\gamma_0$ increases with $p_{\text{select}}$ as expected. Values differ across distributions because $F_v$ differs.

**Note on “mixture”:** at this stage, mixture is treated as a departure affecting the **outcome error** $u$ only (selection error CDF $F_v$ is kept normal). This can be revised if desired.

------------------------------------------------------------------------

## 5. Next implementation steps (pending TA OK)

1)  Fix baseline parameters ($\beta$, $\gamma_{\text{slopes}}$) and a consistent definition of $X$ and $Z$ (including whether to include an exclusion restriction in $Z$).\
2)  Implement the DGP generator and run a smoke test (one scenario, one replication) verifying:
    -   selection rate close to target,
    -   dimensions consistent,
    -   no numerical issues.
3)  Implement estimators + mini-MC sanity runs before scaling to the full grid.

------------------------------------------------------------------------
