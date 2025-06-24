This folder contains the workhorse code to reproduce the analyses in the manuscript

- The helpers/ folder contains two scripts that are sourced in all other files. These scripts contain helper functions.
  - analysis_helpers.R 
  - plot_helpers.R
  
  
- ppe_linear_baseline.R 
  
  - testing a simple baseline yearly model against different seasonal misspecification

- ppe_linear.R 

  - testing the proposed model on simple linear pseudoproxies
  
- ppe_vslite_quantile.R 

  - testing the proposed model on quantile-based VS-Lite pseudoproxies

- ppe_vslite_obs.R 

  - testing the proposed model on observation-based VS-Lite pseudoproxies (i.e. using parameters from Breitenmoser et al 2014)


- ppe_vslite_truncated.R 

  - testing the proposed model on VS-Lite pseudoproxies whose growth is truncated to grow only in DJF

- ppe_vslite_sensitivity.R

  - as in ppe_vslite_truncated, but with multiple experiments to test sensitivies to prior and SNR

- blueoaks_baseline.R

  - a baseline yearly-resolved (DJF) model applied to reconstruct precipitation for the real blue oaks tree ring widths (see also ppe_linear_baseline_eta.R)

- blueoaks_proposed.R 

  - the proposed model fit to the blue oaks data.

- firth_proposed.R 

  - the proposed model fit to the white spruce data.