canonical.index = function (GLM, year, base = BASEYR, year.name = "year", 
                            plots = F, diags = F, line.type = "o", pch = 19, xlim = NULL, 
                            ylim = NULL, rescale = FALSE, add = FALSE, one.line = TRUE, 
                            one.line.lty = 2, one.line.col = "black", xlab = "Fishing year", 
                            ylab = "Index", digits = 4, show.R = FALSE)  {
  res <- as.data.frame(coefficients(GLM))
  index <- regexpr(year.name, row.names(res)) > 0
  X <- res[index, 1]
  V <- summary(GLM)$cov.unscaled[index, , drop = FALSE][, index]
  n <- length(X) + 1
  A <- matrix(-1/n, n, n - 1)
  A[-base, ] <- A[-base, ] + diag(rep(1, n - 1))
  CPUE <- A %*% X
  COV <- sqrt(diag((A %*% V) %*% t(A)))
  index <- exp(CPUE)
  lower.CI <- exp(CPUE - 2 * COV)
  upper.CI <- exp(CPUE + 2 * COV)
  ans <- data.frame(year, index, lower.CI, upper.CI)
  RSQ <- round((GLM$null.deviance - GLM$deviance)/GLM$null.deviance, 
               digits = digits)
  if (is.null(ylim)) {
    ylim <- c(0, ceiling(max(pretty(ans$index))))
  }
  else {
    ylim <- ylim
  }
  if (plots) {
    if (rescale == TRUE) {
      MEAN <- mean(ans$index)
      INDEX <- ans$index/MEAN
      CI.L <- ans$lower.CI/MEAN
      CI.U <- ans$upper.CI/MEAN
    }
    else {
      INDEX <- ans$index
      CI.L <- ans$lower.CI
      CI.U <- ans$upper.CI
    }
    if (!add) {
      plot(ans$year, INDEX, ylim = ylim, xlim = xlim, type = "n", 
           xaxs = "i", yaxs = "i", xlab = xlab, 
           ylab = ylab)
      if (one.line) {
        abline(h = 1, lty = one.line.lty, col = one.line.col)
      }
    }
    points(ans$year, INDEX, pch = pch, type = line.type, 
           bg = one.line.col)
    segments(ans$year, CI.L, ans$year, CI.U)
    if (diags) {
      qqnorm(GLM$residuals)
      qqline(GLM$residuals)
    }
  }
  if (show.R) 
    cat(paste("\n\nOverall r-squared :", RSQ, "\n\n"))
  ans$se <- COV
  ans$cv <- sqrt(exp(COV^2) - 1)
  ans <- list(index = ans, R_squared = RSQ)
  .ans <<- ans
  return(ans)
}
