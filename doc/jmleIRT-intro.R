## ----setup, echo=FALSE, message=FALSE, warning=FALSE--------------------------
  knitr::opts_chunk$set(
    collapse = TRUE,
    comment = "#>",
    fig.width = 10, 
    fig.height = 10
  )
  options(width=80)

  # includes: in_header: "header.html"
  library(jmleIRT)

## ----simulate-data------------------------------------------------------------
N <- 100  # persons
I <- 10   # items
X <- matrix(rbinom(N * I, 1, 0.5), nrow = N)
# randomly set 10% of entries to missing
X[sample(length(X), size = 0.1 * length(X))] <- NA
head(X)


## ----basic-fit----------------------------------------------------------------
fit <- jmle_estimation(
  X, max_iter = 500, conv = 1e-5,
  center = "items", bias_correction = FALSE,
  estimatewle = FALSE, verbose = FALSE
)

str(fit)

## ----bias-correction----------------------------------------------------------
fit_bc <- jmle_estimation(X, center = "items", bias_correction = TRUE)
summary(fit_bc$beta)

## ----bias-compare-------------------------------------------------------------
summary(fit$beta)

## ----wle-example--------------------------------------------------------------
fit_wle <- jmle_estimation(X, center = "items", estimatewle = TRUE)
head(fit_wle$wle_estimate)

## ----missing-data-------------------------------------------------------------
X_miss <- X
X_miss[1:5, 1:2] <- NA
fit_miss <- jmle_estimation(X_miss, center = "items")

length(fit_miss$theta)
length(fit_miss$beta)

