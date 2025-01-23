## code to prepare `Data_Preproc` dataset goes here
## The raw data for the sample can be downloaded from the website:
## https://geodacenter.github.io/data-and-lab/ncovr/.


# —————————————————————————————————————
# Load shapefile data using the st_read function;
# —————————————————————————————————————
NAT <- sf::st_read("/BSTVC/inst/extdata/NAT.shp")

# Since we need to extract a smaller area as our sample map data, we must perform an extraction operation here;
# if the map data you have imported is your study area, then the shapefile data loaded with the st_read function in the previous step is the map data required for modeling.
Florida_Map <- NAT[, c(1:9, 76)] %>% filter(STATE_NAME == "Florida")
# —————————————————————————————————————



# —————————————————————————————————————
# Load the original tabular data, which can be a CSV file or an XLSX file;
# —————————————————————————————————————
# CNTY_NAT <- read.csv("/BSTVC/inst/extdata/NAT.csv")
# CNTY_NAT <- openxlsx::read.xlsx("/BSTVC/inst/extdata/NAT.xlsx")

# since the attribute table of the example shapefile (shp data) we are loading here contains all the variables Y and X, it is also possible to directly load the original tabular data from the shapefile.
# (this method ensures that the order of spatial units in the exported raw tabular data is completely consistent with that of the map data, reducing the likelihood of errors).
CNTY_NAT <- NAT[, 1:69] %>% filter(STATE_NAME == "Florida")
# ————————————————————————————————————-



# —————————————————————————————————————
# Convert the raw data from wide format to long panel data format for BSTVC modeling;
# —————————————————————————————————————
# Since the loaded raw tabular data is in a cross-sectional format (which can be understood as each variable's data for one year being a single column),
# we need to convert it into a spatiotemporal panel data format suitable for modeling (which can be understood as all the data for each variable across all years being a single column).
names(NAT)

# CNTY_NAT[,1:9] represents attribute fields that do not change over time, such as the names and codes of each spatial unit;
# the sample data includes data for four years, so these unchanging attribute fields are duplicated four times, stacked from top to bottom.
Florida_NAT <- rbind.data.frame(CNTY_NAT[, 1:9], CNTY_NAT[, 1:9], CNTY_NAT[, 1:9], CNTY_NAT[, 1:9]) %>%
  dplyr::mutate(Year = rep(c(1960, 1970, 1980, 1990), each = nrow(Florida_Map)))

# Set up an empty data frame to store the  panel data after each transformation
panel <- data.frame()

# The sample data includes a total of 15 variables, so we will perform a loop that iterates 15 times here.
# Each variable's cross-sectional data will be transformed into a spatiotemporal panel format and then combined with the attribute fields that do not change over time.
# 【Core transformation code!!!】
for (i in 1:15) {
  j <- 4 * i + 6
  panel <- as.data.frame(reshape2::melt(
    data = NAT,
    id.vars = colnames(NAT, 9),
    measure.vars = j:j + 3,
    variable.name = paste0("Year", i),
    value.name = stringr::str_extract(names(NAT[j]), "[a-zA-Z]+")[[1]]
  ))
  Florida_NAT <- cbind.data.frame(Florida_NAT, panel[, 10:11])
}

# Remove the redundant year field
names(Florida_NAT)

Florida_NAT <- Florida_NAT[, !names(Florida_NAT) %in% paste0("Year", 1:15)]
# At this point, the transformed Florida_NAT data represents the final spatiotemporal panel data format ready for modeling !!!
# What follows is using the data.check() function to verify whether the order of spatial units in the spatiotemporal panel matches exactly with the order of spatial units in the map data.
# —————————————————————————————————————



# —————————————————————————————————————
# 【Note !!!】
# —————————————————————————————————————
# The variable Y can contain missing values, and a small number of missing values can also exist in variable X. Please represent missing values in the data with NA!!
# Compared to the Geographically Weighted Regression (GWR) and Graphically and Temporally Weighted Regression Model (GTWR) model under the frequentist statistical framework, these are the advantages of our model, namely:

# 1.Our model predicts missing values in variable Y during the modeling process and can also perform spatiotemporal smoothing for all samples of variable Y, which the GTWR model or GWR model cannot do;

# 2.Our model allows for the existence of missing values in variable X, while the GTWR model or GWR model does not.
# —————————————————————————————————————



# —————————————————————————————————————
# Process the count-type Y variable and the binary Y variable in the sample data;
# —————————————————————————————————————
# Users should refer to the information based on their own data.
# Convert count data to integer
Florida_NAT$HC <- round(Florida_NAT$HC)

# Construct binary indicator
Florida_NAT["HW"] <- ifelse(Florida_NAT$HC == 0, 0, 1)
# —————————————————————————————————————
