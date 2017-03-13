#' @include Algorithm.SDM.R
#' @import methods
NULL

setGeneric("data.values", function(obj, Env, na.rm = TRUE) {
  return(standardGeneric("data.values"))
})

setMethod("data.values", "Algorithm.SDM", function(obj, Env, na.rm = TRUE) {
  values <- data.frame(extract(Env, cbind(obj@data$X, obj@data$Y)))

  # Categorical variables as factor
  for (i in seq_len(length(Env@layers))) {
    if (Env[[i]]@data@isfactor) {
      col <- which(names(values) == Env[[i]]@data@names)
      values[, col] <- as.factor(values[, col])
      levels(values[, col]) <- Env[[i]]@data@attributes[[1]]$ID
      if (length(Env[[i]]@data@attributes[[1]]$ID) > 100) {
        warning(paste(names(Env[[i]]), "as more than 100 levels (",
                      length(Env[[i]]@data@attributes[[1]]$ID), ") are you sure to consider it as a factor ?"))
      }
    }
  }

  # Tables binding
  obj@data <- cbind(obj@data, values)

  # NAs removing
  if (na.rm) {
    for (i in seq_len(length(Env@layers))) {
      if (length(which(is.na(obj@data[i + 3]))) > 0) {
        obj@data <- obj@data[-c(which(is.na(obj@data[i + 3]))),
                             ]
      }
    }
  }

  return(obj)
})
