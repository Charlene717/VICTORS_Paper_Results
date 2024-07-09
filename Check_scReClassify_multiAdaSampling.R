
data = sample_sce
label = cellTypes
reducedDimName = "matPCs"
classifier = "svm"
percent = 1
L = 10
prob = FALSE
balance = TRUE
iter = 3


multiAdaSampling <- function(data, label, reducedDimName = NULL, classifier = "svm",
                             percent = 1, L = 10, prob = FALSE, balance = TRUE, iter = 3)
{
  if (is(data, "SingleCellExperiment")) {
    if (is.null(reducedDimName)) {
      stop("reducedDimName parameter cannot be NULL if data is a SingleCellExperiment object")
    }
    mat <- SingleCellExperiment::reducedDim(data, reducedDimName)
  } else {
    mat <- data
  }
  mat <- t(mat)
  if (ncol(mat) != length(label)) {
    stop("Parameter `label` and number of columns of data must be equal")
  }
  if (percent < 0 | percent > 1) {
    stop("Parameter percent must be a value between 0 and 1")
  }
  models <- lapply(seq_len(L), function(l) {
    X <- mat
    Y <- label
    model <- c()
    prob.mat <- c()
    for (i in seq_len(iter)) {
      if (classifier == "rf") {
        model <- randomForest::randomForest(t(X), factor(Y), ntree = 100)
        prob.mat <- stats::predict(model, newdata = t(mat), type = "prob")
      }
      if (classifier == "svm") {
        tmp <- t(X)
        rownames(tmp) <- NULL
        model <- e1071::svm(tmp, factor(Y), probability = TRUE)
        prob.mat <- attr(stats::predict(model, t(mat), decision.values = FALSE, probability = TRUE), "probabilities")
      }
      if (!(classifier %in% c("svm", "rf"))) {
        stop("Classifier ", classifier, " is unknown.")
      }
      X <- c()
      Y <- c()
      XY_update <- lapply(seq_len(ncol(prob.mat)), function(j) {
        voteClass <- prob.mat[label == colnames(prob.mat)[j], ]
        idx <- c()
        if (balance) {
          idx <- sample(seq_len(nrow(voteClass)), size = nrow(voteClass) * percent, replace = TRUE, prob = voteClass[, j])
        } else {
          sampleSize <- round(stats::median(table(label)))
          if (nrow(voteClass) > sampleSize) {
            idx <- sample(seq_len(nrow(voteClass)), size = sampleSize * percent, replace = TRUE, prob = voteClass[, j])
          } else {
            idx <- sample(seq_len(nrow(voteClass)), size = nrow(voteClass) * percent, replace = TRUE, prob = voteClass[, j])
          }
        }
        list(X = mat[, rownames(voteClass)[idx]], Y = label[rownames(voteClass)[idx]])
      })
      cur_X = lapply(XY_update, function(xy) {
        xy$X
      })
      cur_Y = lapply(XY_update, function(xy) {
        xy$Y
      })
      X <- do.call(cbind, cur_X)
      Y <- do.call(c, cur_Y)
    }
    model
  })
  predictMat <- matrix(0, nrow = ncol(mat), ncol = length(table(label)))
  final <- c()
  predmat <- lapply(models, function(m) {
    if (classifier == "svm") {
      tmp <- attr(predict(m, newdata = t(mat), probability = TRUE), "prob")[, names(table(label))]
    } else {
      tmp <- predict(m, newdata = t(mat), type = "prob")[, names(table(label))]
    }
    predictMat <<- predictMat + tmp
  })
  if (prob == TRUE) {
    final <- apply(predictMat, 1, max)
    names(final) <- names(table(label))[apply(predictMat, 1, which.max)]
  } else {
    final <- names(table(label))[apply(predictMat, 1, which.max)]
  }
  return(list(final = final, models = models, prob = predictMat))
}
