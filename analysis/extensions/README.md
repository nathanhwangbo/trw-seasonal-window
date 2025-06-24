This folder contains extensions to the model to facilitate future development.

# Extension 1: Calibration only 

- calibration_only.stan: as in the proposed model, except without any reconstruction components. In particular, the model discrepancy term remains

- calibration_only_noeta.stan: as in calibration_only.stan, except without the model discrepancy components. This is helpful if linearity can be safely assumed

- firth_calib.R: an example application of these models



# Extension 2: Modeling a chronology

- chron_model.stan: as in the proposed model, except the input is a single chronology rather than multiple tree measurements.

- firth_chron.R: an example application of this model, highlighting the dplR package


# Extension 3: Joint preciptiation and temperature modeling

- proposed_model_pt.stan: extending the proposed model to include a joint precipitation and temperature reconstruction

- ppe_linear_pt.R: an idealized pseudoproxy example, using pseudoproxies sensitive to both JJA temperature and DJF preciptiation

- pt_helpers.R: helper functions for ppe_linear_pt.R
