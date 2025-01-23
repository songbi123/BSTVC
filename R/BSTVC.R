#' @title Bayesian Spatiotemporally Varying Coefficients Model
#' @description A class of spatiotemporal regression analysis methods based on Bayesian statistical kernel, whose significant advantage lies in utilizing a "full map" single modeling framework to
#'     uniformly fit the spatiotemporal variations of all local regression coefficients, thereby accurately capturing the spatiotemporal heterogeneity impact of explanatory variables on the target variable, that is, revealing spatiotemporal nonstationarity.
#'
#' @usage BSTVC(
#'     formula = Y ~ ST(X1+X2),
#'     data,
#'     study_map,
#'     Time,
#'     Space,
#'     response_type = c("continuous", "binary", "count"),
#'     threads = 6,
#'     spatial_matrix=NULL,
#'     ...
#'     )
#'
#' @param formula The model formula, with the basic format Y ~ ST(X1 + X2 + X3), is used to specify the target variable and explanatory variables in the modeling. In the BSTVC function, all explanatory variables X must be placed within the ST() symbol, indicating that these explanatory variables have both spatial and temporal nonstationarity (currently only this situation is supported).
#' @param data A data frame containing all variables used in the model. If the data involves multiple time sections, it must be in spatiotemporal panel data format. The order of spatial units under each time section must be completely consistent with the order of geographical units in the shapefile.
#' @param study_map Data in sf format, imported from map data in shapefile (shp) format, which includes unique value fields for each geographical unit as well as geometric attribute information, etc.
#' @param Time A string used to specify the temporal field in the data frame. Only the column name of the temporal field (as a string) needs to be provided.
#' @param Space A string used to specify the unique value field representing spatial units in the data frame. Ensure that the string parameter passed in matches the column name in the data data frame exactly, and that data and study_map both have a field with this string name.
#' @param response_type A string specifying the type of modeling data, currently supporting three application scenarios: "continuous" for continuous Y variables; "binary" for binary Y variables, where 0 indicates non-occurrence and 1 indicates occurrence; "count" for count Y variables.
#' @param threads The number of threads, with a default of 6, which can be set manually.
#' @param spatial_matrix The spatial weight matrix, with a default of 10 nearest neighbors matrix. Users can input their own constructed spatial weight matrix files and support multiple types of spatial weight matrices such as k-nearest, inverse distance, fixed distance, etc.; if the user does not input, this function will help calculate the spatial weight matrix under the QUEEN rule (assigning a value of 1 to units with adjacency relationships, otherwise assigning a value of 0).
#' @param ... NULL
#'
#' @note
#' For detailed guidance on data preprocessing, inspection, model fitting, result output, and result visualization, please visit [GetStart (Chinese)](#) or [GetStart (English)](#).
#'
#'
#' @return
#' The BSTVC function outputs a total of six parts, with specific descriptions for each output section as follows:
#'
#' 1. model.evaluation: Overall evaluation results of the Bayesian model, including commonly used metrics such as deviance information criterion(DIC), watanabe–akaike information criterion(WAIC), and logarithmic score(LS).
#'
#' 2. local.prediction: Local spatiotemporal prediction results for the target variable Y, as well as wide (95%) and narrow (50%) Bayesian credible intervals for each prediction value, used to express uncertainty.
#'
#' 3. summary.random.effects: Random effects results, including spatial random effects and temporal random effects for each explanatory variable, along with their wide (95%) and narrow (50%) Bayesian credible intervals.
#'
#' 4. time.coefficients: Temporal regression coefficients (temporal nonstationarity), in a wide data format, including the temporal regression coefficients of explanatory variables for each time slice and their wide (95%) and narrow (50%) Bayesian credible intervals.
#'
#' 5. space.coefficients: Spatial regression coefficients (spatial nonstationarity), in a wide data format, including the spatial regression coefficients of explanatory variables for each map unit and their wide (95%) and narrow (50%) Bayesian credible intervals.
#'
#' 6. STVPI: Spatiotemporal variance partitioning index (STVPI) calculation results, which quantify the percentage of spatiotemporal contribution of each explanatory variable. This not only provides a ranking of factor importance in terms of spatiotemporal heterogeneity but also assesses relative importance. For more information on this tool, please refer to the literature \href{https://doi.org/10.1016/j.jclepro.2022.133781}{Wan et al., 2022}.
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
#' Bayesian STVC series models: \url{https://chaosong.blog/bayesian-stvc/}
#'
#'
#' @seealso For those interested in INLA and wanting to learn more about Bayesian modeling and inla modeling, you can refer to \code{\link[INLA]{inla}}.
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
#' # Check if the order of sample units for each year in the spatiotemporal panel data table
#' # is consistent with the order of each spatial unit in the map.
#' data <- data.check(
#'   data = Florida_NAT,
#'   study_map = Florida_Map,
#'   Time = "Year",
#'   Space = "FIPS"
#' )
#'
#' # Explanatory variable standardization
#' Florida_NAT[c("DNL","UE","FP","BLK","GI")] <- scale(Florida_NAT[c("DNL","UE","FP","BLK","GI")])
#'
#' model_1 <- BSTVC(HR ~ ST(DNL + UE + FP + BLK + GI),
#'   data = Florida_NAT,
#'   study_map = Florida_Map,
#'   Time = "Year",
#'   Space = "FIPS",
#'   response_type = "continuous",
#'   threads = 6
#' )
#'
#' # class of binary
#' model_2 <- BSTVC(HW ~ ST(DNL + UE + FP + BLK + GI),
#'   data = Florida_NAT,
#'   study_map = Florida_Map,
#'   Time = "Year",
#'   Space = "FIPS",
#'   response_type = "binary",
#'   threads = 6
#' )
#'
#' # class of count
#' model_3 <- BSTVC(HC ~ ST(DNL + UE + FP + BLK + GI),
#'   data = Florida_NAT,
#'   study_map = Florida_Map,
#'   Time = "Year",
#'   Space = "FIPS",
#'   response_type = "count",
#'   threads = 6
#' )
#' }
#'
#'
#' @export




BSTVC <- function(formula = Y ~ ST(X1 + X2), data, study_map, Time, Space, response_type = c("continuous", "binary", "count"), threads = 6, spatial_matrix = NULL, ...) {

  startTime <- Sys.time()

  formula.Y <- strsplit(as.character(formula), " ~ ", fixed = T)[[2]]
  formula.X <- strsplit(sub(".*\\((.*?)\\).*", "\\1", strsplit(as.character(formula), " ~ ", fixed = TRUE)[[3]]), "+", fixed = TRUE)[[1]] %>% trimws()
  num.X <- length(formula.X)
  TimeLength <- length(unique(data[, Time]))
  Spatialareas <- length(data[, Space]) / TimeLength
  data$timeID <- rep(1:TimeLength, each = Spatialareas, times = 1)
  data$spaceID <- rep(1:Spatialareas, each = 1, times = TimeLength)

  data[, paste0("ID.area", 1:num.X)] <- data$spaceID
  data[, paste0("ID.year", 1:num.X)] <- data$timeID


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

  file_path <- system.file("temp_script1.Rc", package = "BSTVC")
  compiler::loadcmp(file_path, envir = rlang::current_env())
  model <- .Z1x3M8pQ(formula.Y,formula.X,data,num.X,link,response_type,adj,threads)

  endTime <- Sys.time()



  DF.cpo <- data.frame(cpo = model$cpo$cpo, pit = model$cpo$pit, failure = model$cpo$failure)

  DF.cpo.best <- subset(DF.cpo, failure == 0, select = c("cpo", "failure"))

  Evaluation_df <- cbind.data.frame(
    ModelName = paste0("BSTVC_Model: Start at ", startTime, ",   end at ", endTime, "."),
    DIC = model$dic$dic,
    eff = model$dic$p.eff,
    LS = -mean(log(DF.cpo.best$cpo)),
    running_time = as.numeric((model$cpu.used / 60)[4]),
    waic = model$waic$waic,
    pd2 = model$waic$p.eff
  )



  if (response_type == "continuous") {
    M.prediction <- cbind(data[Time], data[Space], model$summary.fitted.values) %>%
      `colnames<-`(c(Time, Space, paste0("log(y)_", colnames(model$summary.fitted.values))))

    M.prediction$y <- data[[formula.Y]]
    M.prediction$predict_y <- exp(model$summary.fitted.values$mean)
    M.prediction$fill_y <- ifelse(is.na(M.prediction$y),M.prediction$predict_y,M.prediction$y)
    M.prediction[paste0("predict_", c("0.025quant", "0.25quant", "0.5quant", "0.75quant", "0.975quant"))] <-
      exp(M.prediction[, paste0("log(y)_", c("0.025quant", "0.25quant", "0.5quant", "0.75quant", "0.975quant"))])
    M.prediction <- M.prediction[, !colnames(M.prediction) %in% "log(y)_mode"]

  } else if (response_type == "binary") {
    M.prediction <- cbind(data[Time], data[Space], model$summary.fitted.values) %>%
      `colnames<-`(c(Time, Space, paste0("y_", colnames(model$summary.fitted.values))))

    M.prediction$y <- data[[formula.Y]]
    M.prediction$predict_y <- model$summary.fitted.values$mean
    M.prediction <- M.prediction[, !colnames(M.prediction) %in% "y_mode"]

  } else if (response_type == "count") {
    M.prediction <- cbind(data[Time], data[Space], model$summary.fitted.values) %>%
      `colnames<-`(c(Time, Space, paste0("y_", colnames(model$summary.fitted.values))))

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
      paste0("SD.for.ID.area", 1:num.X),
      paste0("SD.for.ID.year", 1:num.X)
    ), ]
    summary.random.SD["random effects"] <- c(
      rownames(summary.random.SD)[1],
      paste0("SD.space.", formula.X[1:num.X]),
      paste0("SD.time.", formula.X[1:num.X])
    )
  } else if (response_type == "binary") {
    summary.random.SD <- summary.random.SD[c(
      paste0("SD.for.ID.area", 1:num.X),
      paste0("SD.for.ID.year", 1:num.X)
    ), ]
    summary.random.SD["random effects"] <- c(
      paste0("SD.space.", formula.X[1:num.X]),
      paste0("SD.time.", formula.X[1:num.X])
    )
  } else if (response_type == "count") {
    summary.random.SD <- summary.random.SD[c(
      paste0("SD.for.ID.area", 1:num.X),
      paste0("SD.for.ID.year", 1:num.X)
    ), ]
    summary.random.SD["random effects"] <- c(
      paste0("SD.space.", formula.X[1:num.X]),
      paste0("SD.time.", formula.X[1:num.X])
    )
  } else {
    summary.random.SD <- summary.random.SD[c(
      rownames(summary.random.SD)[1],
      paste0("SD.for.ID.area", 1:num.X),
      paste0("SD.for.ID.year", 1:num.X)
    ), ]
    summary.random.SD["random effects"] <- c(
      rownames(summary.random.SD)[1],
      paste0("SD.space.", formula.X[1:num.X]),
      paste0("SD.time.", formula.X[1:num.X])
    )
  }

  rownames(summary.random.SD) <- NULL




  df.model.random.time <- do.call(rbind.data.frame, model$summary.random[paste0("ID.year", 1:num.X)]) %>%
    dplyr::mutate(
      time_index = rep(1:TimeLength, each = 1, times = num.X),
      explain_variable = rep(formula.X, each = TimeLength, times = 1)
    ) %>%
    `rownames<-`(NULL)

  df.model.random.time <- df.model.random.time[, !(colnames(df.model.random.time) %in% c("ID", "mode", "kld"))]





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
      paste0("Precision.for.ID.area", 1:num.X),
      paste0("Precision.for.ID.year", 1:num.X)
    )]
  } else if (response_type == "binary") {
    samp.vars.main.random <- samp.vars[, c(
      paste0("Precision.for.ID.area", 1:num.X),
      paste0("Precision.for.ID.year", 1:num.X)
    )]
  } else if (response_type == "count") {
    samp.vars.main.random <- samp.vars[, c(
      paste0("Precision.for.ID.area", 1:num.X),
      paste0("Precision.for.ID.year", 1:num.X)
    )]
  } else {
    samp.vars.main.random <- samp.vars[, c(
      paste0("Precision.for.the.", response_type, ".observations"),
      paste0("Precision.for.ID.area", 1:num.X),
      paste0("Precision.for.ID.year", 1:num.X)
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

  samp.vars.main.random[paste0("time.", formula.X[1:num.X])] <-
    samp.vars.main.random[paste0("Precision.for.ID.year", 1:num.X)] /
      samp.vars.main.random$model.original.total.SD



  samp.vars.main.random$space.time.all <-
    rowSums(samp.vars.main.random[, c(
      paste0("Precision.for.ID.area", 1:num.X),
      paste0("Precision.for.ID.year", 1:num.X)
    )]) / samp.vars.main.random$model.original.total.SD

  samp.vars.main.random$space.all <-
    rowSums(samp.vars.main.random[, paste0("Precision.for.ID.area", 1:num.X)]) / samp.vars.main.random$model.original.total.SD

  samp.vars.main.random$time.all <-
    rowSums(samp.vars.main.random[, paste0("Precision.for.ID.year", 1:num.X)]) / samp.vars.main.random$model.original.total.SD


  if (num.X > 1) {
    for (k in 1:num.X) {
      samp.vars.main.random[paste0("space.time.", formula.X[k])] <-
        rowSums(samp.vars.main.random[, c(
          paste0("Precision.for.ID.area", k),
          paste0("Precision.for.ID.year", k)
        )]) /
          samp.vars.main.random$model.original.total.SD
    }
  } else if (num.X == 1) {
    samp.vars.main.random[paste0("space.time.", formula.X[1])] <-
      rowSums(samp.vars.main.random[, c(
        "Precision.for.ID.area1",
        "Precision.for.ID.year1"
      )]) /
        samp.vars.main.random$model.original.total.SD
  }


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
    random.LGM.contribution.95CI <- random.LGM.contribution.95CI[-c(1:(2 * num.X + 2)), ]
  } else if (response_type == "binary") {
    random.LGM.contribution.95CI <- random.LGM.contribution.95CI[-c(1:(2 * num.X + 1)), ]
  } else if (response_type == "count") {
    random.LGM.contribution.95CI <- random.LGM.contribution.95CI[-c(1:(2 * num.X + 1)), ]
  } else {
    random.LGM.contribution.95CI <- random.LGM.contribution.95CI[-c(1:(2 * num.X + 1)), ]
  }




  return.list <- list(
    model.evaluation = Evaluation_df,
    local.prediction = M.prediction,
    summary.random.effects = summary.random.SD,
    time.coefficients = df.model.random.time,
    space.coefficients = df.model.random.space,
    STVPI = random.LGM.contribution.95CI
  )
  return(return.list)
}
