## Source code for Learning and Predicting from Dynamic Models for COVID-19 Patient Monitoring"
Repository for paper "Learning and Predicting from Dynamic Models for COVID-19 Patient Monitoring"

* *00_loadlibs.R*: setup file directories and load required pacakges 

* Code folder

   * generatePrediction
   
      1. *multinom_results_alldata.R*: run multinomial regression models for prospective and retrospective appraoch using all available data, and generate coefficient tables.
      
      2. *jtgaussian_results_alldata.R*: run multinomial regression models for prospective and retrospective appraoch, and generate Figure 1.
      
      3. *model_perform_forward.R*: run overall and cross-validated prospective linear mixed effects model; run cross-validated prospective multinomial model; generate prospective method predictions for each patient provided patients are hospitalized for 0, 2, 4, 8 days (baseline days).
      
      4. *model_perform_backward.R*: run overall and cross-validated retrospective linear mixed effects model; run cross-validated retrospective multinomial model; generate retrospective method predictions for each patient provided patients are hospitalized for 0, 2, 4, 8 days.
      
   * modelChecking
   
      1. *mixed_model_checking.R*: produce residuals versus fitted deciles plot and residuals Q-Q plot against standard Gaussian for prospective and retrospective linear mixed models.
      
      2. *discrim_auc_script.R*: read predictions from *model_perform_forward.R* and *model_perform_backward.R* and calculate AUC from predictions for each baseline day.
      
      3. *calib_chisqstats_script.R*: read prediction outputs and calculate chisquare statistics for calibration power.
      
      4. *produce_calib_discrim_figures.R*: read outputs from *discrim_auc_script.R* and *calib_chisqstats_script.R* and generate Figure 4.
