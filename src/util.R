#' normalize the log probabilities
#' 
#' @param logw a vector that contains the unnormalized log probabilities
#' @return w a vector that contains the normalized probabilities
#'  
normalize.logw <- function(logw){
  # guard against underflow or overflow
  c <- max(logw)
  w <- exp(logw-c)
  # normalize the probabilities
  w <- w / sum(w)
  return(w)
}

#' convert null results in the mat file to a data frame
#' 
#' @param null.path a string that indicates the path of mat file
#' @return null.df a data frame that contains the results in the mat file
#'  
null.mat2df <- function(null.path){
  null.results <- R.matlab::readMat(null.path)
  
  null.df <- data.frame(
    theta0 <- c(null.results$theta0.vec),
    h <- c(null.results$h.vec),
    logw.step1 <- c(null.results$logw.step1.vec),
    logw.step2 <- c(null.results$logw.step2.vec)
  )
  names(null.df) <- c("theta0","h","logw.step1","logw.step2")
  
  # normalize the log importance weights
  null.df$posp.step1 <- normalize.logw(null.df$logw.step1)
  null.df$posp.step2 <- normalize.logw(null.df$logw.step2)
  
  return(null.df)
}

#' convert gsea results in the mat file to a data frame
#' 
#' @param gsea.path a string that indicates the path of mat file 
#' @return gsea.df a data frame that contains the results in the mat file 
#' 
gsea.mat2df <- function(gsea.path){
  gsea.data <- R.matlab::readMat(gsea.path)
  
  gsea.df <- data.frame(
    id <- c(gsea.data$path.id),
    name <- unlist(gsea.data$path.name),
    source <- unlist(gsea.data$path.sr),
    database <- unlist(gsea.data$path.base),
    numgene <- c(gsea.data$path.numg),
    numsnps <- c(gsea.data$path.nums),
    log10.bf <- log10(c(gsea.data$path.bf)),
    theta.mean <- gsea.data$theta.est[, 1],
    theta.95lb <- gsea.data$theta.est[, 3],
    theta.95ub <- gsea.data$theta.est[, 4]
  )
  names(gsea.df) <- c("id","name","source","database","numgene","numsnps",
                      "log10.bf","theta.mean","theta.95lb","theta.95ub")
  
  # convert the names of pathways to ASCII
  gsea.df$name <- as.character(gsea.df$name)
  gsea.df$name <- iconv(gsea.df$name, from = "latin1", to = "UTF-8")
  
  return(gsea.df)
}