# approxstats

## Overview
Software companion for the methodology proposed in *"Data-driven stabilizations of goodness-of-fit tests"* (Fernández-de-Marcos and García-Portugués, 2021).

## Usage

### $(n, alpha)$-stabilization

In order to fit a $(n, \alpha)$-model to the ratios $T_{\infty; \alpha}/T_{\n; \alpha}$, use the script [fit_quantiles_nalpha.R](). Execute all steps, specfying the statistic you want to fit. 

**IMPORTANT!** Pre-computed Monte Carlo quantiles must be available at the folder */distributions/*.

```R
# Specify a statistic from avail_cir_tests
avail_cir_tests
statistic <- "KS"
```

When fitting a model, you have several input options. Check [src/fit_T_Tn_stable_ratios.R]() to see a detailed description. However, if you want to reproduce the results from the paper, you might be interested in the following:

```R
T_Tn_fitted <- fit_T_Tn(. . .,
                        fit.direction = c("both", "forward"),
                        weights.fun = weight.fun.2,
                        lambda = 2, mu = 2)
```

- **fit.direction**: Specifies the optimal-search direction, can be "all", "both", or "forward".
- **weights.fun**: Weight function that can be chosen from [src/weights.R](), numbered in the same order as in *Appendix A*.
- **lambda**: Maximum power for $n$-like predictors, of the form (n^{-lambda/2}). Minimum required: 2. Check *Appendix B* for more detailed information.
- **mu**: Maximum power for $\alpha$-like predictors, of the form (n^{-mu/2}). Minimum required: 2. Check *Appendix B* for more detailed information.

### $(n, \alpha, p)$-stabilization

## Data application in astronomy
The data application in Fernández-de-Marcos and García-Portugués (2021) can be reproduced through the script []().

## References
Fernández-de-Marcos, A., and García-Portugués, E. (2021). Data-driven stabilizations of goodness-of-fit tests. *arXiv:2112.01808*. https://arxiv.org/abs/2112.01808.