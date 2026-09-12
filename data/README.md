# Data inputs

This directory contains the two preprocessed Human Mortality Database input objects used by the analysis:

```text
HMD_Dx.RData
HMD_Ex.RData
```

`HMD_Dx.RData` contains death-count inputs and `HMD_Ex.RData` contains exposure inputs. Population selection and analysis restrictions are applied by the R code using the metadata and settings in `config/`.

The underlying source data are from the Human Mortality Database (HMD), https://www.mortality.org/. Users should consult HMD for documentation, citation guidance and terms of use.
