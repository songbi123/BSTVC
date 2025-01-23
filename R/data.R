#'
#' Murder data by county in the U.S.
#'
#' @description
#' Homicide and selected socioeconomic characteristics of counties in the continental United States. Data for four decennial census years: 1960, 1970, 1980 and 1990.
#'     This dataset is commonly used to analyze trends in violent crime, particularly murder rates, in different geographic areas and to explore the relationship between social, economic, and environmental factors and crime rates. We have processed this raw data as follows:
#'
#' 1. Considering the size of the dataset and the fact that this data is used solely as sample data, we have chosen to extract only the county data for Florida from the original dataset to serve as our example. This selection is aimed at ensuring the manageability of the dataset while providing a specific geographic focus for our analysis.
#'
#' 2. Furthermore, since the "Homicide count" field in the dataset is calculated as an average over three years centered on the current census year, it includes decimal values. To align this field with the characteristics of a count variable, we have rounded the values to the nearest whole number.
#'
#' 3. Additionally, to meet the model's requirement for binary (0 or 1) data, indicating the occurrence or non-occurrence of homicides, we have computed a new "HW" field based on the processed "Homicide count" field.
#'
#'
#' @format A data frame with 268 rows and 26 variables:
#' \describe{
#'  \item{NAME}{County name}
#'  \item{STATE_NAME}{State name}
#'  \item{STATE_FIPS}{State fips code (character)}
#'  \item{CNTY_FIPS}{County fips code (character)}
#'  \item{FIPS}{Combined state and county fips code (character)}
#'  \item{STFIPS}{State fips code (numeric)}
#'  \item{COFIPS}{County fips code (numeric)}
#'  \item{FIPSNO}{Fips code as numeric variable (State fips*1000+County fips)}
#'  \item{SOUTH}{Dummy variable for Southern counties (South = 1)}
#'
#'  \item{Year}{Data year, that is, 1960, 1970, 1980, 1990.}
#'  \item{HR}{Homicide Rate per 100,000 (numerator is a 3 year average centered on the decennial census year, e.g., 1959, 1960, and 1961). Individual deaths are aggregated to the county level according to the decedent's county of residence.}
#'  \item{HC}{Homicide count (3 year average centered on decennial census year, e.g., 1959, 1960, and 1961). Individual deaths are aggregated to the county level according to the decedent's county of residence.}
#'  \item{HW}{Binary data (0, 1), indicating whether a homicide has occurred. If HC is 0, this field is 0; if HC is not 0, this field is 1.}
#'  \item{PO}{Population of each county for the decennial census year}
#'  \item{RD}{Resource Deprivation/Affluence Component (principal component composed of percent black, log of median family income, gini index of family income inequality, percent of families female headed (percent of families single parent for 1960), and percent of families below poverty (percent of families below $3,000 for 1960) (See Land et al., 1990))}
#'  \item{PS}{Population Structure Component (principal component composed of the log of population and the log of population density) (See Land et al., 1990)}
#'  \item{UE}{Percent of labor is unemployed}
#'  \item{DV}{Percent of males 14 and over who are divorced (aged 15 and over for 1980 and 1990)}
#'  \item{MA}{Median age}
#'  \item{POL}{Population logged}
#'  \item{DNL}{Population density logged}
#'  \item{MFIL}{Median family income logged}
#'  \item{FP}{Percent of families below poverty (percent of families below $3,000 for 1960)}
#'  \item{BLK}{Percent black}
#'  \item{GI}{Gini index of family income inequality}
#'  \item{FH}{Percent of families female headed (percent of families single parent for 1960)}
#' }
#'
#' @source S. Messner, L. Anselin, D. Hawkins, G. Deane, S. Tolnay, R. Baller (2000). An Atlas of the Spatial Patterning of County-Level Homicide, 1960-1990. Pittsburgh, PA,National Consortium on Violence Research (NCOVR).
#' Data can be downloaded from the \href{https://geodacenter.github.io/data-and-lab/ncovr/}{GeoDa Data and Lab}.
#'
#' @references Baller, R., L. Anselin, S. Messner, G. Deane and D. Hawkins (2001). Structural covariates of US county homicide rates: incorporating spatial effects. Criminology 39, 561-590.
#'
"Florida_NAT"




#'
#' Map of County Districts in Florida, USA.
#'
#' The shapefile data for the 67 counties in Florida, United States, which includes not only unique value fields and county names as attribute information for each county, but also geometric information for each county.
#'
#' @format A object of class \code{\link[sf]{st_read}} with 67 rows and 10 variables:
#' \describe{
#'  \item{NAME}{County name}
#'  \item{STATE_NAME}{State name}
#'  \item{STATE_FIPS}{State fips code (character)}
#'  \item{CNTY_FIPS}{County fips code (character)}
#'  \item{FIPS}{Combined state and county fips code (character)}
#'  \item{STFIPS}{State fips code (numeric)}
#'  \item{COFIPS}{County fips code (numeric)}
#'  \item{FIPSNO}{Fips code as numeric variable (State fips*1000+County fips)}
#'  \item{SOUTH}{Dummy variable for Southern counties (South = 1)}
#'  \item{geometry}{Geometric information for each county}
#' }
#'
#' @source S. Messner, L. Anselin, D. Hawkins, G. Deane, S. Tolnay, R. Baller (2000). An Atlas of the Spatial Patterning of County-Level Homicide, 1960-1990. Pittsburgh, PA,National Consortium on Violence Research (NCOVR).
#' Data can be downloaded from the \href{https://geodacenter.github.io/data-and-lab/ncovr/}{GeoDa Data and Lab}.
#'
#' @references Baller, R., L. Anselin, S. Messner, G. Deane and D. Hawkins (2001). Structural covariates of US county homicide rates: incorporating spatial effects. Criminology 39, 561-590.
#'
"Florida_Map"
