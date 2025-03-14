#' @title Data Checking of Modeling
#'
#' @description
#' Based on the key field check, whether the order of spatial units in the spatiotemporal panel data is completely consistent with the order of geographical units in the map. If there is at least one year with an inconsistent order, this function helps to re-match the data.
#'
#' @usage
#' data.check(data,study_map,Time = "",Space = "")
#'
#' @param data A data frame containing all variables used in the model (including a key field representing unique spatial units, and the type of the key field must be consistent in both data and study_map), with n years of data in a spatiotemporal panel format.
#' @param study_map Data in sf format, imported from shapefile (shp) map data. It includes unique value fields for each geographical unit as well as geometric attribute information, etc.
#' @param Time A string used to specify the temporal field in the data frame. Ensure that the string parameter passed in matches the column name in the data frame exactly.
#' @param Space A string used to specify the unique value field representing spatial units in the data frame. Make sure that the string parameter passed in matches the column name in the data frame exactly, and that both data and study_map have a field with this string name.
#'
#' @return If the order of spatial units in the input data for each year is completely consistent with the order of geographical units in the regional map, only a prompt will be output, meaning users can proceed to the next step of modeling based on the input data. If there is an inconsistency, a newly re-matched dataset will be output.
#'
#' @note Note that the number of geographical units in the map data (e.g., 30 districts) should be completely consistent with the number of units in the data; the unique value fields for spatial units (Space) in data and study_map must be of the same format, both numeric, character types, etc.
#'
#' @seealso For methods on converting cross-sectional data (data for each year as a column) into spatiotemporal panel data, refer to the \code{\link[reshape2]{melt}} function in the reshape2 package or the \code{\link[tidyr]{pivot_longer}} function in the tidyr package. For precautions during the data checking process, please refer to the getstart guide in the vignettes folder.
#'
#' @importFrom dplyr arrange sym
#'
#' @export
#'
#' @examples
#' data(Florida_NAT)
#' data(Florida_Map)
#'
#' # Note: pay attention to check the attributes of each field.
#' str(Florida_NAT)
#' str(Florida_Map)
#' # Florida_Map$FIPS <- as.numeric(Florida_Map$FIPS)
#'
#' newdata <- data.check(Florida_NAT,
#'   study_map = Florida_Map,
#'   Time = "Year",
#'   Space = "FIPS"
#' )


data.check <- function(data, study_map, Time, Space) {
  unique.years <- unique(data[[Time]])
  num.years <- length(unique.years)
  studymap <- as.data.frame(study_map)

  check.year <- vector()
  for (Y in 1:num.years) {
    result <- all(data[sum(1,nrow(studymap)*(Y-1)):(nrow(studymap)*Y),Space] == studymap[,Space])
    check.year <- c(check.year, result)
  }

  if (!all(check.year)) {
    cat("Error: The number of spatial units per year in the data does not match the number of geographical units in the study_map!\n")

    cat("\n...At least one unique value field of a temporal node does not match; currently in the process of data matching...\n")

    data.neworder <- dplyr::arrange(data, !!sym(Time), !!sym(Space))
    data <- data.neworder

    data$Code <- data[[Space]]
    study_map$CODE <- study_map[[Space]]

    order_InfirstCircle <- match(study_map$CODE, data$Code)

    order_matrix <- matrix(nrow = nrow(study_map), ncol = num.years)

    for (i in 1:num.years) {
      order_matrix[, i] <- as.vector(order_InfirstCircle + nrow(study_map) * (i - 1))
    }

    for (j in 1:num.years) {
      if (j == 1) {
        orderfinal = order_matrix[,j]
      }
      else {
        orderfinal=rbind(orderfinal, order_matrix[,j])
      }
    }

    orderfinal_Transpose <- t(orderfinal)

    orderfinal.vector <- as.vector(orderfinal_Transpose)

    Orderdf <- data[orderfinal.vector, ]

    Orderdf$timeID <- rep(1:num.years, each = nrow(study_map), times = 1)

    Orderdf$spaceID <- rep(1:nrow(study_map), each = 1, times = num.years)

    check.org <- vector()
    studymap <- as.data.frame(study_map)

    for (years in unique.years) {
      check <- Orderdf[Orderdf[Time] == years, "Code"] == studymap["CODE"]
      check.org <- c(check.org, check)
    }


    if (all(check.org)) {
      data <- Orderdf
      cat("\nSuccess! The data has been matched. Please use the new data for subsequent modeling!")
      return(data)
    }
    else{
      cat("Error: The data did not match successfully. Please check the data or the types and formats of the fields!")
    }
  }
  else if(all(check.year)){
    cat("The unique value field of the spatial unit completely matches, and you can proceed to the next step of modeling!")
    return(data)
  }
}
