# TIVAP Fibrin Sheath Prediction

This repository contains the source code and web application for an interpretable machine-learning model developed to predict fibrin sheath formation after totally implantable venous access port (TIVAP) implantation.

## Overview

The model integrates catheter spatial measurements with routinely available haematological and coagulation-related variables to provide individualized risk estimates of TIVAP-related fibrin sheath formation.

The simplified model includes six predictors:

* Distance from the catheter tip to the carina
* Distance from the catheter tip to the clavicle
* Distance from the catheter tip to the right atrium
* Activated partial thromboplastin time
* Mean platelet volume
* Monocyte count

The final prediction model was developed using an XGBoost linear booster and evaluated in training, internal testing, and temporal validation cohorts.

## Repository contents

* `app.R`: Source code for the Shiny web application
* `model/`: Trained model object and supporting files
* `R/`: Functions used for preprocessing and prediction
* `example/`: Example input data
* `README.md`: Project documentation
* `LICENSE`: Open-source license

## Web application

The interactive prediction tool is available at:

https://portsheathrisk.shinyapps.io/FS_predict/

## Requirements

The application was developed in R version 4.4.2.

Required R packages may include:

```r
shiny
xgboost
dplyr
ggplot2
```

Install the required packages using:

```r
install.packages(c(
  "shiny",
  "xgboost",
  "dplyr",
  "ggplot2"
))
```

## Running the application locally

Clone or download this repository, open the project directory in RStudio, and run:

```r
shiny::runApp()
```

## Data availability

The patient-level data used to develop and validate the model are not publicly available because of patient privacy and institutional restrictions. De-identified data may be available from the corresponding author upon reasonable request and subject to approval by the institutional ethics committee.

No identifiable patient information is included in this repository.

## Intended use

This model is intended for research purposes and exploratory clinical risk assessment. It should not replace clinical judgment or established institutional protocols. Further multicentre prospective validation is required before routine clinical implementation.

## Citation

If you use this code or web application, please cite:

Chen C, Wu C, Zong Y, Bao Y, Yang X, Zhao S, Gu W. Spatial–Biological Interaction Drives Fibrin Sheath Formation: An Interpretable Machine Learning Study in Totally Implantable Venous Access Ports. Manuscript submitted for publication.

The citation information will be updated after publication.

## Contact

Corresponding author:

Weiwei Gu
Department of Interventional and Vascular Surgery
Affiliated Hospital of Nantong University
Nantong, China
Email: [5201189@ntu.edu.cn](mailto:5201189@ntu.edu.cn)

## License

This project is licensed under the MIT License.
