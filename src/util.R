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
  
  # extract log 10 BFs of gene set enrichment
  # the following order ensures that log(bf), if exisiting, will be used
  if ("path.bf" %in% names(gsea.data)) {
    log10.bf <- log10(c(gsea.data$path.bf))
  }
  if ("path.lbf" %in% names(gsea.data)) {
    log10.bf <- c(gsea.data$path.lbf) / log(10) 
  }
  
  # create output data frame
  gsea.df <- data.frame(
    id <- c(gsea.data$path.id),
    name <- unlist(gsea.data$path.name),
    source <- unlist(gsea.data$path.sr),
    database <- unlist(gsea.data$path.base),
    numgene <- c(gsea.data$path.numg),
    numsnps <- c(gsea.data$path.nums),
    log10.bf <- log10.bf,
    theta.mean <- gsea.data$theta.est[, 1],
    theta.95lb <- gsea.data$theta.est[, 3],
    theta.95ub <- gsea.data$theta.est[, 4]
  )
  names(gsea.df) <- c("id","name","source","database","numgene","numsnps",
                      "log10.bf","theta.mean","theta.95lb","theta.95ub")
  
  # convert the names of pathways to ASCII
  gsea.df$name <- as.character(gsea.df$name)
  gsea.df$name <- iconv(gsea.df$name, from = "latin1", to = "UTF-8")
  
  # remove gene sets whose source and database are 'NA'
  keep.index <- (gsea.df$database!="NA") & (gsea.df$source!="NA")
  gsea.df <- gsea.df[keep.index, ]
  
  # sort rows by BF values, from largest to smallest
  gsea.df <- plyr::arrange(gsea.df, -log10.bf)
  
  return(gsea.df)
}

gtex.mat2df <- function(gtex.path){
  gtex.data <- R.matlab::readMat(gtex.path)
  
  # extract log 10 BFs of gene set enrichment
  if ("path.bf" %in% names(gtex.data)) {
    log10.bf <- log10(c(gtex.data$path.bf))
  }
  if ("path.lbf" %in% names(gtex.data)) {
    log10.bf <- c(gtex.data$path.lbf) / log(10) 
  }
  
  # create output data frame
  gtex.df <- data.frame(
    id <- c(gtex.data$path.id),
    numsnps <- c(gtex.data$path.nums),
    log10.bf <- log10.bf,
    theta.mean <- gtex.data$theta.est[, 1],
    theta.95lb <- gtex.data$theta.est[, 3],
    theta.95ub <- gtex.data$theta.est[, 4]
  )
  names(gtex.df) <- c("id","numsnps","log10.bf","theta.mean","theta.95lb","theta.95ub")
  
  return(gtex.df)
}

gtex.mat2df.round2 <- function(gtex.path){
  gtex.data <- R.matlab::readMat(gtex.path)
  
  # extract log 10 BFs of gene set enrichment
  log10.bf <- c(gtex.data$path.lbf) / log(10)
  
  # create gene set id
  if ("path.id" %in% names(gtex.data)) {
    path.id <- as.character(c(gtex.data$path.id))
  }
  
  if ("path.ts" %in% names(gtex.data)) {
    path.ts <- unlist(gtex.data$path.ts)
    path.id <- gsub("^.*?_","",path.ts) 
  }
  
  # create output data frame
  gtex.df <- data.frame(
    id <- path.id,
    numsnps <- c(gtex.data$path.nums),
    log10.bf <- log10.bf,
    theta0.mean <- gtex.data$theta0.est[, 1],
    theta0.95lb <- gtex.data$theta0.est[, 3],
    theta0.95ub <- gtex.data$theta0.est[, 4],
    theta.mean <- gtex.data$theta.est[, 1],
    theta.95lb <- gtex.data$theta.est[, 3],
    theta.95ub <- gtex.data$theta.est[, 4]
  )
  names(gtex.df) <- c("id","numsnps","log10.bf",
                      "theta0.mean","theta0.95lb","theta0.95ub",
                      "theta.mean","theta.95lb","theta.95ub")
  
  # combine mean with 95 C.I.
  printed.theta0.mean <- sprintf("%.3f", round(gtex.df$theta0.mean,digits=3))
  printed.theta0.95lb <- sprintf("%.3f", round(gtex.df$theta0.95lb,digits=3))
  printed.theta0.95ub <- sprintf("%.3f", round(gtex.df$theta0.95ub,digits=3))
  
  printed.theta.mean <- sprintf("%.3f", round(gtex.df$theta.mean,digits=3))
  printed.theta.95lb <- sprintf("%.3f", round(gtex.df$theta.95lb,digits=3))
  printed.theta.95ub <- sprintf("%.3f", round(gtex.df$theta.95ub,digits=3))
  
  gtex.df$theta0 <- paste0(printed.theta0.mean,", [",printed.theta0.95lb,", ",printed.theta0.95ub,"]")
  gtex.df$theta <- paste0(printed.theta.mean,", [",printed.theta.95lb,", ",printed.theta.95ub,"]")
  
  gtex.df.output <- gtex.df[, c("id","numsnps","log10.bf","theta0","theta")]
  
  # sort rows by BF values, from largest to smallest
  gtex.df.output <- plyr::arrange(gtex.df.output, -log10.bf)
  
  return(gtex.df.output)
  
}

gsea.mat2df.round2 <- function(gsea.path){
  gsea.data <- R.matlab::readMat(gsea.path)
  
  # extract log 10 BFs of gene set enrichment
  if ("path.bf" %in% names(gsea.data)) {
    log10.bf <- log10(c(gsea.data$path.bf))
  }
  if ("path.lbf" %in% names(gsea.data)) {
    log10.bf <- c(gsea.data$path.lbf) / log(10) 
  }
  
  # create output data frame
  gsea.df <- data.frame(
    id <- c(gsea.data$path.id),
    name <- unlist(gsea.data$path.name),
    source <- unlist(gsea.data$path.sr),
    database <- unlist(gsea.data$path.base),
    numgene <- c(gsea.data$path.numg),
    numsnps <- c(gsea.data$path.nums),
    log10.bf <- log10.bf,
    theta0.mean <- gsea.data$theta0.est[, 1],
    theta0.95lb <- gsea.data$theta0.est[, 3],
    theta0.95ub <- gsea.data$theta0.est[, 4],
    theta.mean <- gsea.data$theta.est[, 1],
    theta.95lb <- gsea.data$theta.est[, 3],
    theta.95ub <- gsea.data$theta.est[, 4]
  )
  names(gsea.df) <- c("id","name","source","database","numgene","numsnps",
                      "log10.bf","theta0.mean","theta0.95lb","theta0.95ub",
                      "theta.mean","theta.95lb","theta.95ub")
  
  # convert the names of pathways to ASCII
  gsea.df$name <- as.character(gsea.df$name)
  gsea.df$name <- iconv(gsea.df$name, from = "latin1", to = "UTF-8")
  
  # combine mean with 95 C.I.
  printed.theta0.mean <- sprintf("%.3f", round(gsea.df$theta0.mean,digits=3))
  printed.theta0.95lb <- sprintf("%.3f", round(gsea.df$theta0.95lb,digits=3))
  printed.theta0.95ub <- sprintf("%.3f", round(gsea.df$theta0.95ub,digits=3))
  
  printed.theta.mean <- sprintf("%.3f", round(gsea.df$theta.mean,digits=3))
  printed.theta.95lb <- sprintf("%.3f", round(gsea.df$theta.95lb,digits=3))
  printed.theta.95ub <- sprintf("%.3f", round(gsea.df$theta.95ub,digits=3))
  
  gsea.df$theta0 <- paste0(printed.theta0.mean,", [",printed.theta0.95lb,", ",printed.theta0.95ub,"]")
  gsea.df$theta <- paste0(printed.theta.mean,", [",printed.theta.95lb,", ",printed.theta.95ub,"]")
  
  gsea.df.output <- gsea.df[, c("id","name","source","database","numgene","numsnps","log10.bf","theta0","theta")]
  
  # remove gene sets whose source and database are 'NA'
  keep.index <- (gsea.df.output$database!="NA") & (gsea.df.output$source!="NA")
  gsea.df.output <- gsea.df.output[keep.index, ]
  
  # sort rows by BF values, from largest to smallest
  gsea.df.output <- plyr::arrange(gsea.df.output, -log10.bf)
  
  return(gsea.df.output)
}

