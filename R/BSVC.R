#' @title Bayesian Spatially Varying Coefficients Model
#' @description The Bayesian Spatially Varying Coefficients (BSVC) model is a spatially simplified version of the BSTVC model,
#'     designed specifically for identifying variable relationships with spatial heterogeneity, also known as spatial nonstationarity.
#'     Its advantage lies in integrating the "full map" single modeling framework of BSTVC, ensuring that the fitted local spatial regression coefficients are directly comparable.
#'
#' @usage BSVC(
#'     formula =Y~S(X1+X2),
#'     data,
#'     study_map,
#'     Space,
#'     response_type = c("continuous", "binary", "count"),
#'     threads = 6,
#'     spatial_matrix=NULL,
#'     ...
#'     )
#'
#'
#'
#' @param formula The model formula, with the basic format Y ~ S(X1 + X2 + X3), is used to specify the target variable and explanatory variables in the modeling. In the BSVC function, all explanatory variables X must be placed within the S() symbol, indicating that the explanatory variables X exhibit spatial nonstationarity.
#' @param data A data frame containing all variables used in the model. This data frame should only include information that varies in the spatial dimension, i.e., it should represent a time slice of spatiotemporal panel data, such as data for a specific year or specific month. The order of spatial units in each year in the data must be completely consistent with the order of geographical units in the shapefile.
#' @param study_map Data in sf format, imported from map data in shapefile (shp) format, which includes unique value fields for each geographical unit as well as geometric attribute information, etc.
#' @param Space A string used to specify the unique value field representing spatial units in the data frame. Ensure that the string parameter passed in matches the column name in the data data frame exactly, and that data and study_map both have a field with this string name.
#' @param response_type A string specifying the type of modeling data, currently supporting three application scenarios: "continuous" for continuous Y variables; "binary" for binary Y variables, where 0 indicates non-occurrence and 1 indicates occurrence; "count" for count Y variables.
#' @param threads The number of threads, with a default of 6, which can be set manually.
#' @param spatial_matrix  The spatial weight matrix, with a default of a 10-neighbor matrix. Users can input their own constructed spatial weight matrix files, and support multiple types of spatial weight matrices such as k-nearest, inverse distance, fixed distance, etc.; if the user does not input, this function will help calculate the binary spatial weight matrix under the QUEEN rule (assigning a value of 1 to units with adjacency relationships, otherwise assigning a value of 0).
#' @param ... NULL
#'
#' @note
#' For detailed guidance on data preprocessing, inspection, model fitting, result output, and result visualization, please visit [GetStart (Chinese)](#) or [GetStart (English)](#).
#'
#'
#' @return
#' The BSVC function outputs a total of five parts, with specific descriptions for each output section as follows:
#'
#' 1. model.evaluation: Overall evaluation results of the Bayesian model, including commonly used metrics such as deviance information criterion(DIC), watanabe–akaike information criterion(WAIC), and logarithmic score(LS).
#'
#' 2. local.prediction: Local prediction results for the target variable Y, along with wide (95%) and narrow (50%) Bayesian credible intervals for each predicted value, used to express uncertainty.
#'
#' 3. summary.random.effects: Random effects results, including spatial random effects for each explanatory variable, as well as wide (95%) and narrow (50%) Bayesian credible intervals.
#'
#' 4. space.coefficients: Spatial regression coefficients (SCs), in a wide data format, including spatial regression coefficients for each map unit and their wide (95%) and narrow (50%) Bayesian credible intervals, used to characterize spatial nonstationarity.
#'
#' 5. STVPI: Spatiotemporal variance partitioning index (STVPI) calculation results, which quantify the spatiotemporal contribution percentage of each explanatory variable. This will output the spatial contribution percentage for each explanatory variable X. For more information on this tool, please refer to the literature \href{https://doi.org/10.1016/j.jclepro.2022.133781}{Wan et al., 2022}.
#'
#'
#' @author \href{https://chaosong.blog/}{Chao Song}; Xianteng Tang;
#'
#' @references
#' Chao Song, et al., Exploring spatiotemporal nonstationary effects of climate factors on hand, foot, and mouth disease using Bayesian Spatiotemporally Varying Coefficients (STVC) model in Sichuan, China, Sci. Total Environ. 648 (2019) 550–560.
#'
#' Chao Song, et al., Spatiotemporally Varying Coefficients (STVC) model: a Bayesian local regression to detect spatial and temporal nonstationarity in variables relationships, ANN GIS 26 (3) (2020) 277–291.
#'
#' Chao Song, et al., (2022). Spatiotemporal disparities in regional public risk perception of COVID-19 using Bayesian Spatiotemporally Varying Coefficients (STVC) series models across Chinese cities. International Journal of Disaster Risk Reduction, 103078.
#'
#' Wan Qin, Tang Zhangying, Pan Jay, Xie Mingyu, Wang Shaobin, Yin Hao, Li Junmin, Liu Xin, Yang Yang, & Song Chao. (2022). Spatiotemporal heterogeneity in associations of national population ageing with socioeconomic and environmental factors at the global scale. Journal of Cleaner Production, 373, 133781.
#'
#' Bayesian STVC series models: https://chaosong.blog/bayesian-STVC/.
#'
#'
#' @seealso For those interested in INLA and wanting to learn more about Bayesian modeling and inla modeling, you can refer to \code{\link[INLA]{inla}}.
#'
#'
#'
#'
#' @importFrom spdep poly2nb nb2mat
#' @importFrom stats as.formula quantile
#' @importFrom INLA inla inla.hyperpar.sample
#' @importFrom dplyr mutate %>%
#' @importFrom stringr str_extract_all
#' @importFrom scales percent_format
#' @importFrom compiler loadcmp
#' @importFrom rlang current_env
#'
#' @examples \donttest{
#' # For detailed data processing procedures, modeling steps, and precautions,
#' # please refer to the document at "vignettes/getstart.Rmd".
#' data(Florida_Map)
#' data(Florida_NAT)
#'
#' # Core field inspection, such as the "Time" field, "Space" field
#' str(Florida_Map)
#' str(Florida_NAT)
#'
#' # Extract a temporal slice from the data to serve as the input data
#' # for the spatially varying coefficient model.
#' data <- Florida_NAT[Florida_NAT$Year == 1990, ]
#'
#' # Check if the order of sample units for each year in the spatiotemporal panel data table
#' # is consistent with the order of each spatial unit in the map.
#' newdata <- data.check(
#'   data = data,
#'   study_map = Florida_Map,
#'   Time = "Year",
#'   Space = "FIPS"
#' )
#'
#' # "Explanatory variable standardization"
#' newdata[c("DNL", "UE", "FP", "BLK", "GI")] <- scale(newdata[c("DNL", "UE", "FP", "BLK", "GI")])
#'
#' model_1 <- BSVC(HR ~ S(DNL + UE + FP + BLK + GI),
#'   data = newdata,
#'   study_map = Florida_Map,
#'   Space = "FIPS",
#'   response_type = "continuous"
#' )
#'
#' # class of binary
#' model_2 <- BSVC(HW ~ S(DNL + UE + FP + BLK + GI),
#'   data = newdata,
#'   study_map = Florida_Map,
#'   Space = "FIPS",
#'   response_type = "binary"
#' )
#'
#' # class of count
#' model_3 <- BSVC(HC ~ S(DNL + UE + FP + BLK + GI),
#'   data = newdata,
#'   study_map = Florida_Map,
#'   Space = "FIPS",
#'   response_type = "count"
#' )
#' }
#'
#'
#' @export





BSVC <- function(formula = Y ~ S(X1 + X2), data, study_map, Space, response_type = c("continuous", "binary", "count"), threads = 6, spatial_matrix = NULL, ...) {

  startTime <- Sys.time()

  formula.Y <- strsplit(as.character(formula), " ~ ", fixed = T)[[2]]
  formula.X <- strsplit(sub(".*\\((.*?)\\).*", "\\1", strsplit(as.character(formula), " ~ ", fixed = TRUE)[[3]]), "+", fixed = TRUE)[[1]] %>% trimws()
  num.X <- length(formula.X)
  Spatialareas <- length(data[, Space])
  data$spaceID <- rep(1:Spatialareas, each = 1, times = 1)
  data[, paste0("ID.area", 1:num.X)] <- data$spaceID


  if (is.null(spatial_matrix)) {
    NeighboursList <- spdep::poly2nb(study_map, queen = T)
    adj <- spdep::nb2mat(NeighboursList, style = "B", zero.policy = TRUE)
  } else {
    adj <- spatial_matrix
  }


  TaregtDataColumn <- data[, formula.Y]
  data$TargetY <- TaregtDataColumn
  sum(is.na(TaregtDataColumn))
  link <- rep(NA, length(TaregtDataColumn))
  link[which(is.na(TaregtDataColumn))] <- 1


  file_path <- system.file("temp_script2.Rc", package = "BSTVC")
  compiler::loadcmp(file_path, envir = rlang::current_env())
  model <- .Jv_9gA7xU4(num.X,formula.X,data,formula.Y,link,response_type,adj,threads)

  endTime <- Sys.time()



  DF.cpo <- data.frame(cpo = model$cpo$cpo, pit = model$cpo$pit, failure = model$cpo$failure)

  DF.cpo.best <- subset(DF.cpo, failure == 0, select = c("cpo", "failure"))

  Evaluation_df <- cbind.data.frame(
    ModelName = paste0("BSVC_Model: Start at ", startTime, ",   end at ", endTime, "."),
    DIC = model$dic$dic,
    eff = model$dic$p.eff,
    LS = -mean(log(DF.cpo.best$cpo)),
    running_time = as.numeric((model$cpu.used / 60)[4]),
    waic = model$waic$waic,
    pd2 = model$waic$p.eff
  )




  if (response_type == "continuous") {
    M.prediction <- cbind(data[Space], model$summary.fitted.values) %>%
      `colnames<-`(c(Space, paste0("log(y)_", colnames(model$summary.fitted.values))))

    M.prediction$y <- data[[formula.Y]]
    M.prediction$predict_y <- exp(model$summary.fitted.values$mean)
    M.prediction$fill_y <- ifelse(is.na(M.prediction$y),M.prediction$predict_y,M.prediction$y)
    M.prediction[paste0("predict_", c("0.025quant", "0.25quant", "0.5quant", "0.75quant", "0.975quant"))] <-
      exp(M.prediction[, paste0("log(y)_", c("0.025quant", "0.25quant", "0.5quant", "0.75quant", "0.975quant"))])
    M.prediction <- M.prediction[, !colnames(M.prediction) %in% "log(y)_mode"]

  } else if (response_type == "binary") {
    M.prediction <- cbind(data[Space], model$summary.fitted.values) %>%
      `colnames<-`(c(Space, paste0("y_", colnames(model$summary.fitted.values))))

    M.prediction$y <- data[[formula.Y]]
    M.prediction$predict_y <- model$summary.fitted.values$mean
    M.prediction <- M.prediction[, !colnames(M.prediction) %in% "y_mode"]

  } else if (response_type == "count") {
    M.prediction <- cbind(data[Space], model$summary.fitted.values) %>%
      `colnames<-`(c(Space, paste0("y_", colnames(model$summary.fitted.values))))

    M.prediction$y <- data[[formula.Y]]
    M.prediction$predict_y <- model$summary.fitted.values$mean
    M.prediction <- M.prediction[, !colnames(M.prediction) %in% "y_mode"]
  } else {
    stop("Error: currently does not support other types.")
  }




  n.sample <- 100000

  summary.random.SD.samp <- data.frame(sqrt(1 / INLA::inla.hyperpar.sample(n.sample, model)))
  summary.random.SD <- as.data.frame(t(apply(
    summary.random.SD.samp,
    2,
    function(x) {
      stats::quantile(x, c(0.025, 0.25, 0.5, 0.75, 0.975))
    }
  )))

  rownames(summary.random.SD) <- gsub("Precision", "SD", rownames(summary.random.SD))

  if (response_type == "continuous") {
    summary.random.SD <- summary.random.SD[c(
      rownames(summary.random.SD)[1],
      paste0("SD.for.ID.area", 1:num.X)
    ), ]
    summary.random.SD["random effects"] <- c(
      rownames(summary.random.SD)[1],
      paste0("SD.space.", formula.X[1:num.X])
    )
  } else if (response_type == "binary") {
    summary.random.SD <- summary.random.SD[paste0("SD.for.ID.area", 1:num.X), ]
    summary.random.SD["random effects"] <- paste0("SD.space.", formula.X[1:num.X])
  } else if (response_type == "count") {
    summary.random.SD <- summary.random.SD[paste0("SD.for.ID.area", 1:num.X), ]
    summary.random.SD["random effects"] <- paste0("SD.space.", formula.X[1:num.X])
  } else {
    summary.random.SD <- summary.random.SD[c(
      rownames(summary.random.SD)[1],
      paste0("SD.for.ID.area", 1:num.X)
    ), ]
    summary.random.SD["random effects"] <- c(
      rownames(summary.random.SD)[1],
      paste0("SD.space.", formula.X[1:num.X])
    )
  }

  rownames(summary.random.SD) <- NULL




  Mapindex <- unique(data[[Space]])

  df.model.random.space <- do.call(cbind.data.frame, model$summary.random[paste0("ID.area", 1:num.X)]) %>%
    `colnames<-`(paste(rep(formula.X, each = 10, times = 1),
      rep(colnames(model$summary.random$ID.area1), each = 1, times = num.X),
      sep = "_"
    ))
  df.model.random.space[names(data[Space])] <- Mapindex

  character <- paste0(formula.X, rep(c("_mode", "_kld"), each = num.X))

  df.model.random.space <- df.model.random.space[, !(colnames(df.model.random.space) %in% character)]





  samp.vars <- data.frame(sqrt(1 / INLA::inla.hyperpar.sample(n.sample, model)))

  if (response_type == "continuous") {
    samp.vars.main.random <- samp.vars[, c(
      "Precision.for.the.lognormal.observations",
      paste0("Precision.for.ID.area", 1:num.X)
    )]
  } else if (response_type == "binary") {
    samp.vars.main.random <- samp.vars[, paste0("Precision.for.ID.area", 1:num.X)]
  } else if (response_type == "count") {
    samp.vars.main.random <- samp.vars[, paste0("Precision.for.ID.area", 1:num.X)]
  } else {
    samp.vars.main.random <- samp.vars[, c(
      paste0("Precision.for.the.", response_type, ".observations"),
      paste0("Precision.for.ID.area", 1:num.X)
    )]
  }



  samp.vars.main.random$model.original.total.SD <- rowSums(samp.vars.main.random)

  if (response_type == "continuous") {
    samp.vars.main.random$residual <-
      samp.vars.main.random$Precision.for.the.lognormal.observations /
        samp.vars.main.random$model.original.total.SD
  } else if (response_type == "binary") {
    samp.vars.main.random
  } else if (response_type == "count") {
    samp.vars.main.random
  } else {
    samp.vars.main.random$residual <-
      samp.vars.main.random[paste0("Precision.for.the.", response_type, ".observations")] /
        samp.vars.main.random$model.original.total.SD
  }


  samp.vars.main.random[paste0("space.", formula.X[1:num.X])] <-
    samp.vars.main.random[paste0("Precision.for.ID.area", 1:num.X)] /
      samp.vars.main.random$model.original.total.SD


  samp.vars.main.random$space.all <-
    rowSums(samp.vars.main.random[, paste0("Precision.for.ID.area", 1:num.X)]) / samp.vars.main.random$model.original.total.SD



  random.LGM.contribution.95CI <- as.data.frame(t(apply(
    samp.vars.main.random,
    2,
    function(x) {
      stats::quantile(x, c(0.025, 0.25, 0.5, 0.75, 0.975))
    }
  ))) %>%
    `colnames<-`(paste0("STVPI_", colnames(.)))

  random.LGM.contribution.95CI["random effects"] <- rownames(random.LGM.contribution.95CI)
  colnames(random.LGM.contribution.95CI) <- c("STVPI_2.5%", "STVPI_25%", "STVPI_mean", "STVPI_75%", "STVPI_97.5%", "random effects")

  rownames(random.LGM.contribution.95CI) <- NULL
  if (response_type == "continuous") {
    random.LGM.contribution.95CI <- random.LGM.contribution.95CI[-c(1:(num.X + 2)), ]
  } else if (response_type == "binary") {
    random.LGM.contribution.95CI <- random.LGM.contribution.95CI[-c(1:(num.X + 1)), ]
  } else if (response_type == "count") {
    random.LGM.contribution.95CI <- random.LGM.contribution.95CI[-c(1:(num.X + 1)), ]
  } else {
    random.LGM.contribution.95CI <- random.LGM.contribution.95CI[-c(1:(num.X + 1)), ]
  }



  return.list <- list(
    model.evaluation = Evaluation_df,
    local.prediction = M.prediction,
    summary.random.effects = summary.random.SD,
    space.coefficients = df.model.random.space,
    STVPI = random.LGM.contribution.95CI
  )
  return(return.list)
}
