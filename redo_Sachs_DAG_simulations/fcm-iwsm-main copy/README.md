# Falsifying causal models via conditional independence testing

This repository contains the code for reproducing the results in [1].

# Reproducibility

The results can be reproduced by running `make all`, or executing the following
steps:
```bash
cd code 
Rscript --vanilla dependencies.R # Installs dependencies
Rscript --vanilla graphs.R # Produces Figure 1
Rscript --vanilla main.R # Produces Table 1
```

# References

[1] L. Kook (2025). Falsifying causal models via conditional independence
testing. In: Proceedings of the 39th International Workshop on Statistical
Modelling (IWSM), Limerick, Ireland.
