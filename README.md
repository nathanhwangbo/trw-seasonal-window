Code to reproduce analysis from "Simultaneous optimal target season estimation and local climate reconstruction using tree-ring widths"

To get a sense of the modeling workflow, I recommend starting with /analysis/r_code/ppe_linear.R

- pre-processed time series used in the study are stored in /data/

To better understand how the model works on the `stan` side of things, I think /analysis/stan_code/baseline_model is a good place to start (to get a sense of how a yearly resolved model works), before moving to the proposed model which estimates a seasonal window.

Questions or comments can be directed to nhwangbo@ucla.edu


