# iSEXS_cellReports2019
This is the code necessary for the downloading, preprocessing, and analysis within the publication Bongen et al. 2019 Cell Reports


### How to use this repo
  1. Source **generalGEMfunx.R** for preprocessing and general purpose functions
  2. Source **quickPlots.R** for ggplot2 wrapper functions that use a MetaIntegrator Dataset object as input
  3. Use the appropriate scripts within **downloadData** to create the Metaintegrator Dataset objects you need
  4. To run the meta-analysis of Discovery or Validation cohorts, use **metaAnalysis.R**
  5. To generate the plots and analyses seen in the paper, use **analysisFigures.R**
  

