# approxstats

## Overview
Software companion for the methodology proposed in *"Data-driven stabilizations of goodness-of-fit tests"* (Fernández-de-Marcos and García-Portugués, 2021).

## Data

Data can be found in [RunMyCode.org]().

- `distributions.zip` contains precomputed quantiles (1e7 Monte Carlo samples) for statistics $D_n$, $W^2_n$, $A^2_n$, $P^{AD}\_{n; p}$, $P^{CvM}\_{n; p}$, and $N_{n; p}$; and hyperspherical $n=500$ approximation of the asymptotic quantiles in `/asymptotic`. This data is necessary to fit a $(n, \alpha, d)$-stabilization model.

## Fit a $(n, \alpha, d)$-stabilization model for a statistic $T_n$

### $(n, \alpha)$-stabilization

In order to fit a $(n, \alpha)$-model to the ratios $T_{\infty; \alpha}/T_{n; \alpha}$, use the script [fit_quantiles_nalpha.R](https://github.com/afernandezdemarcos/approxstats/blob/main/fit_quantiles_nalpha.R). Execute all steps, specfying the statistic you want to fit. 

**IMPORTANT!** Pre-computed Monte Carlo quantiles must be available at the folder */distributions/*.

```R
# Specify a statistic from avail_cir_tests
avail_cir_tests
statistic <- "KS"
```

When fitting a model, you have several input options. Check [src/fit_T_Tn_stable_ratios.R](https://github.com/afernandezdemarcos/approxstats/blob/main/src/fit_T_Tn_stable_ratios.R) to see a detailed description. However, if you want to reproduce the results from the paper, you might be interested in the following:

```R
T_Tn_fitted <- fit_T_Tn(. . .,
                        fit.direction = c("both", "forward"),
                        weights.fun = weight.fun.2,
                        lambda = 2, mu = 2)
```

- **fit.direction**: Specifies the optimal-search direction, can be "all", "both", or "forward".
- **weights.fun**: Weight function that can be chosen from [src/weights.R](https://github.com/afernandezdemarcos/approxstats/blob/main/src/weights.R), numbered in the same order as in *Appendix A*.
- **lambda**: Maximum power for $n$-like predictors, of the form $n^{-\lambda/2}$. Minimum required: 2. Check *Appendix B* for more detailed information.
- **mu**: Maximum power for $\alpha$-like predictors, of the form $\alpha^{-\mu/2}$. Minimum required: 2. Check *Appendix B* for more detailed information.

### $(n, \alpha, p)$-stabilization

In order to fit a $(n, \alpha, p)$-model to the ratios $T_{\infty; \alpha; p}/T_{n; \alpha; p}$, use the script [fit_quantiles_nalphap.R](https://github.com/afernandezdemarcos/approxstats/blob/main/fit_quantiles_nalphap.R). Execute all steps, specfying the statistic you want to fit. 

**IMPORTANT!** Asymptotic quantiles are approximated by $n=500$ quantiles, which must be available at the folder */distributions/asymptotic/*.

```R
# Specify a statistic from avail_sph_tests
avail_sph_tests
statistic <- "PCvM"
```

It works similar to the $(n, \alpha)$-stabilization, except for the data preparation step, which is made directly in the main script. In addition, the **model.reference** is set to the specific form presented in (8).

## Approximation assessment

In order to assess the performance of the fit compared to Monte Carlo and particular approximation methods for each statistic, both in terms of accuracy and efficiency, the following scripts are to be used.

**Accuracy**
- 

**Computation efficiency**
- [exec_time/exec_time_comparison.R](https://github.com/afernandezdemarcos/approxstats/blob/main/exec_time/exec_time_comparison.R): Saves and analyzes execution times for Algorithm 1, Monte Carlo and particular approximation methods in `/results`.

## Uniformity tests

In order to perform circular and hyperspherical tests using the $(n, \alpha, p)$ stabilization and Algorithm 1 p-value approximation, use `unif_test_mod` from [src/unif_test.R](https://github.com/afernandezdemarcos/approxstats/blob/main/src/unif_test.R).

```R
n <- 10
samp_cir <- r_unif_cir(n = n)

V_n <- unif_test_mod(Theta = samp_cir, 
                          statistic = "Kuiper")
W2_n <- unif_test_mod(Theta = samp_cir, 
                          statistic = "Watson")

samp_sph <- r_unif_sph(n = n, p = 3)
PCvM_n <- unif_test_mod(Theta = samp_sph, 
                          statistic = "PCvM")
PAD_n <- unif_test_mod(Theta = samp_sph, 
                          statistic = "PAD")
Bakshaev_n <- unif_test_mod(Theta = samp_sph, 
                          statistic = "Bakshaev")
```

## Data application in astronomy

The data application can be reproduced through the following scripts. `rotasym` package contains sunspots data.

**Uniformity test computation**
- [sunspots.R](https://github.com/afernandezdemarcos/approxstats/blob/main/sunspots/sunspots.R): Computes the (n, alpha)-stabilized p-value approximation, and saves the results and execution times.
- [sunspots_MC.R](https://github.com/afernandezdemarcos/approxstats/blob/main/sunspots/sunspots_MC.R): Computes the Monte Carlo (5e3 samples) p-value approximation, and saves its execution times.

**Analysis**

Once the uniformity tests are computed, the analysis can be performed in:
- [analysis_sunspots.R](https://github.com/afernandezdemarcos/approxstats/blob/main/sunspots/analysis_sunspots.R): It builds the whole analysis presented at Figure 4.
- [exec_time_sunspots.R](https://github.com/afernandezdemarcos/approxstats/blob/main/sunspots/exec_time_sunspots.R): Comparison of execution times between Monte Carlo and (n, alpha)-stabilization.

The results are stored in `/sunspots/results`.


## References

Fernández-de-Marcos, A., and García-Portugués, E. (2021). Data-driven stabilizations of goodness-of-fit tests. *arXiv:2112.01808*. https://arxiv.org/abs/2112.01808.