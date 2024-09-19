# Biosignature Detection with Py GC-MS Data using Machine Learning

## The R-scripts are created for the paper

### Detecting Biosignatures in Complex Molecular Mixtures from pyrolysis Gas – Chromatography – Mass Spectrometry Data using Machine Learning

Grethe Hystad<sup>1</sup>, H. James Cleaves II <sup>2,3,4</sup>, Collin A. Garmon<sup>5</sup>, Michael L. Wong<sup>6,7</sup>, Anirudh Prabhu<sup>6</sup>, George D. Cody<sup>6</sup>, and Robert M. Hazen<sup>6</sup>

 *1. Department of Mathematics and Statistics, Purdue University Northwest, Hammond, IN, 46323, USA.*
 
 *2. Department of Chemistry, Howard University, Washington, D.C. 20059, USA.*
 
 *3. Earth Life Science Institute, Tokyo Institute of Technology, Tokyo, Japan.*
 
 *4. Blue Marble Space Institute for Science, Seattle, WA 98104, USA.*
 
 *5. Department of Mathematical Sciences, Purdue University Fort Wayne, Fort Wayne, IN, 46805, USA. *
 
 *6. Earth and Planets Laboratory, Carnegie Science, Washington, DC 20015, USA.*
 
 *7. NHFP Sagan Fellow, NASA Hubble Fellowship Program, Space Telescope Science Institute, Baltimore, MD 21218, USA.*

## Introduction
Three-dimensional (scan number /mass-to-charge ratio/intensity) data from biotic and abiotic samples are obtained by pyrolysis-gas chromatography mass spectrometry. The R-scripts created are for preprocessing these data and to use machine learning to predict whether a sample is biotic or abiotic. Nested resampling is used to obtain an estimate for the prediction performance of the model. The pattern of features that are "important" in distinguishing the abiotic from the biotic species are then determined and shown graphically.
The following machine learning classification methods are used: random forest, logistic regression with elastic net penalty, support vector machines (SVM), and eXtreme Gradient Boosting (XGBoost). The Benjamini-Hochberg procedure is used for multiple hypothesis testing. 

## Data
The first 134 pyr-GC-MS data given in "Cleavesetal.pyrGCMSData.zip", can be found at https://accounts.osf.io/login?service=https://osf.io/embh8/ with reference:

Cleaves, H. J. (2023). A robust, agnostic molecular biosignature based on machine learning (Version 1) [Dataset]. OSF. DOI 10.17605/OSF.IO/EMBH8

The next 16 pyr-GC-MS data are given in "NEWOSFData.zip" (currently set to private).
The R Markdown files can be found in RPubs at: https://rpubs.com/ghystad

## Licence
The application is released under GNU GPL version 3 license.

## Author of the R-scripts
Grethe Hystad

## Sessioninfo

R version 4.3.3 (2024-02-29 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 11 x64 (build 22621)

Matrix products: default


locale:
[1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8    LC_MONETARY=English_United States.utf8
[4] LC_NUMERIC=C                           LC_TIME=English_United States.utf8    

time zone: America/Chicago
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices datasets  utils     methods   base     

loaded via a namespace (and not attached):
 [1] compiler_4.3.3    fastmap_1.1.1     cli_3.6.2         htmltools_0.5.8.1 tools_4.3.3       yaml_2.3.8        rmarkdown_2.26   
 [8] knitr_1.46        xfun_0.43         digest_0.6.35     rlang_1.1.3       renv_1.0.7        evaluate_0.23    

 Cite as: [![DOI](https://zenodo.org/badge/860000284.svg)](https://zenodo.org/doi/10.5281/zenodo.13799340)
 

 
