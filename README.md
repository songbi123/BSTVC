
<!-- README.md is generated from README.Rmd. Please edit that file -->

# BSTVC <img src="man/Figure/R_logo.png" alt="BSTVC" align="right" height="160"/>

<!-- badges: start -->

![CRAN](https://www.r-pkg.org/badges/version/BSTVC)
![Github](https://img.shields.io/badge/publish-2025_01_28-edddab)
![download2](https://cranlogs.r-pkg.org/badges/grand-total/BSTVC)
![download](https://cranlogs.r-pkg.org/badges/BSTVC)
![GitHub](https://img.shields.io/github/license/songbi123/BSTVC)
[![DOI](https://zenodo.org/badge/DOI/10.1016/j.ijdrr.2022.103078.svg)](https://doi.org/10.1016/j.ijdrr.2022.103078)

<!-- badges: end -->

**Spatiotemporal heterogeneous perspective for analyzing influencing
factors, identifying key drivers, and making dynamic predictions, all
within a unified ‘full-map’ framework.**

<!-- The BSTVC package offers a comprehensive and unified "full-map" geographic modeling framework designed to accurately capture spatiotemporal disparities in variable relationships. Its primary goal is to uncover spatiotemporal heterogeneous impacts of multiple explanatory variables on the target variable, i.e., spatiotemporal nonstationarity (Song et al., 2019, 2020, 2022; Wan et al., 2022).  -->
<!-- Our BSTVC package is user-friendly, catering to the in-depth needs of professionals while lowering the barriers to complex Bayesian modeling. This makes advanced Bayesian local spatiotemporal regression methods accessible to a broader user community, enabling easier analysis and interpretation of complex spatiotemporal panel data. It is applicable across a wide range of disciplines, including but not limited to public health, medical geography, environmental health, health economics, and social medicine (Song and Tang, 2025). -->

## Installation

**- Install the `BSTVC` R package**

The package is currently in the internal testing phase. At present, it
only supports local installation from GitHub.

``` r
# Install using the devtools package
# install.packages("devtools")
devtools::install_github("songbi123/BSTVC")

# Install using the remotes package
# install.packages("remotes")
remotes::install_github("songbi123/BSTVC")
```

**- Install the dependency package - the `INLA` R package**

When installing the `BSTVC` package in RStudio, the system will prompt
you to install additional R packages that come with it. However, since
`INLA` is a larger package, installing the `BSTVC` package might lead to
a failure. To avoid this issue, we provide a separate method for
installing the `INLA` package for your reference.

If the installation of the `BSTVC` package in the previous step failed,
please install the `BSTVC` package after successfully installing the
`INLA` package. If you have successfully installed the `INLA` package
while installing the `BSTVC` package, you can skip this step.

``` r
## To install the INLA package, more information can be found at <https://www.r-inla.org/download-install>.
# Extend the overtime duration to 5 minutes
options(timeout = 300)
install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)

## Another approach, using an older version of INLA, is to download the compressed package of the INLA R package to your local machine and then proceed with the installation.
```

## Features & Advantages

The `BSTVC` R package is designed to provide a comprehensive suite of
functionalities for advanced spatiotemporal heterogeneous analysis.
Here’s what our package can do for you:

- **Targeting multiple types of response variables**: It supports three
  mainstream types of response variables: continuous (log-Gaussian
  regression), binary (logistic regression), and count (Poisson
  regression), accommodating various analytical scenarios.

- **Detecting spatiotemporal heterogeneous impact mechanisms**: By
  fitting spatiotemporal regression coefficients, it reveals local
  spatiotemporal differences between explanatory variables (X) and
  response variables (Y), facilitating an in-depth analysis of
  context-specific patterns and exploring the impact mechanisms brought
  by spatiotemporal heterogeneity.

- **Identifying spatiotemporal driving factors**: On the basis of
  identifying spatiotemporal heterogeneous impact mechanisms, it
  clarifies key driving factors by calculating the spatiotemporal
  explainable percentage, providing strong evidence for geographical
  spatiotemporal attribution.

- **Improving spatiotemporal prediction accuracy**: Considering the
  spatiotemporal heterogeneity of local variable relationships, it
  significantly improves model fitting and prediction accuracy, which
  can be used for spatiotemporal missing value imputation,
  spatiotemporal smoothing, and future forecasting.

- **Bayesian model assessment**: It provides a comprehensive evaluation
  of Bayesian regression models, including model fitting (DIC, WAIC),
  complexity (pd), and prediction accuracy (LS) indicators, helping
  users fully understand model performance.

- **Rich visualization outputs**: It provides a variety of
  spatiotemporal visualization tools and codes to help users intuitively
  understand model results, enhance the interpretability of data
  analysis, and promote innovation in your applied research.

Bayesian STVC model is a powerful analytical tool with many advantages
that other similar tools lack, such as **a *“full-map” modeling
framework*, *parameter uncertainty*, *friendliness to missing values*,
and *support for more spatial weight matrices***, among others.

## Usage Guide

To help you quickly and fully get started with our R package for complex
data analysis, we have prepared several detailed and comprehensive usage
guides, as follows:

| Guide | Details |
|----|----|
| **User’s Guide for the BSTVC R Package** | This usage guide covers detailed example operations and important considerations for each key step, including data import, inspection, preprocessing, model fitting, result output and result visualization. You can view it in the `GetStart.Rmd` document under the `vignettes` folder, but it’s in R markdown format. <br><br>If you want to download the help document in PDF format, please click [here](https://github.com/songbi123/BSTVC/raw/songbi123-useguides/GetStart.pdf), the filename is [GetStart-English.pdf](https://github.com/songbi123/BSTVC/raw/songbi123-useguides/GetStart.pdf). At the same time, to meet the needs of Chinese users, we have also provided a Chinese version of the usage guide, which can be downloaded and saved locally by visiting [用户手册-中文版.pdf](https://github.com/songbi123/BSTVC/raw/songbi123-useguides/GetStart-Chinese.pdf). |
| **Modeling Data Processing Guide** | This usage guide demonstrates how to import the types of data required for the model and how to transform the raw data into the spatiotemporal panel data format that can be processed by the BSTVC model. The R code for achieving data processing for modeling is located in the `Data_Preproc.R` file under the `data-raw` folder. |

In the near future, we will continue to refine our documentation and
provide new help documents.

## Changelog

View detailed changelog: [CHANGELOG.md](./CHANGELOG.md)

## Contact

We welcome and encourage user contributions, including reporting issues,
requesting new features, or submitting code changes. If you encounter
any problems when using the BSTVC package or need further assistance,
you can get support through the following means:

1.  **GitHub issues**: Report issues or request new features in
    theGitHub repository, please visit
    [Issues](https://github.com/songbi123/BSTVC/issues).
2.  **Email contact**: <tangxxxxt@163.com>(Tang Xianteng, related to R
    package usage); <chaosong.gis@gmail.com> (Song Chao, related to
    statistical theory)
3.  **Bayesian STVC model**: <https://chaosong.blog/bayesian-stvc/>

Copyright **©HEOA-West China Health and Medical Geography Research
Group**

If you are a WeChat user, you are welcome to scan the QR code to follow
our research group’s official account: **HealthGeography**

<img src="man/Figure/wechat.png" alt="WeChat" align="center" height="180"/>

## Reference

- **\[Bayesian STVC series models\]** Song, Chao, Yin, Hao, Shi, Xun,
  Xie, Mingyu, Yang, Shujuan, Zhou, Junmin, Wang, Xiuli, Tang,
  Zhangying, Yang, Yili, & Pan, Jay. (2022). Spatiotemporal disparities
  in regional public risk perception of COVID-19 using Bayesian
  Spatiotemporally Varying Coefficients (STVC) series models across
  Chinese cities. *International Journal of Disaster Risk Reduction*,
  77, 103078.

- **\[STVPI\]** Wan, Qin, Tang, Zhangying, Pan, Jay, Xie, Mingyu, Wang,
  Shaobin, Yin, Hao, Li, Junmin, Liu, Xin, Yang, Yang, & Song, Chao.
  (2022). Spatiotemporal heterogeneity in associations of national
  population ageing with socioeconomic and environmental factors at the
  global scale. *Journal of Cleaner Production*, 373, 133781.

- Song, Chao, Shi, Xun, & Wang, Jinfeng. (2020). Spatiotemporally
  Varying Coefficients (STVC) model: a Bayesian local regression to
  detect spatial and temporal nonstationarity in variables
  relationships. *Annals of GIS*, 26(3), 277-291.

- Song, Chao, Shi, Xun, Bo, Yanchen, Wang, Jinfeng, Wang, Yong, & Huang,
  Dacang. (2019). Exploring Spatiotemporal Nonstationary Effects of
  Climate Factors on Hand, Foot, and Mouth Disease Using Bayesian
  Spatiotemporally Varying Coefficients (STVC) Model in Sichuan, China.
  *Science of The Total Environment*, 648, 550-560.
