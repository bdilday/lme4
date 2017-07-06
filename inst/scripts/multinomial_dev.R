
#' @export
test_multinomial_data <- function(nlim=NULL) {
  BProDRA::generate_model_df(nlim=nlim, year=2016)  
}

#' @export
test_binomial_data <- function(nlim=NULL) {
  model_df <-  BProDRA::generate_model_df(nlim=nlim, year=2016)  
  cc <- which(model_df$ev$outcome == 1)
  model_df$ev[-cc,]$outcome <- 0
  model_df
  
}

test_nomial_model <- function(model_df, family_name='binomial', ...) {
  mod1 <- lme4::glmer(outcome ~ (1|bid) + (1|pid) + (1|sid), 
                      data=model_df$ev, 
                      nAGQ = 0, 
                      family=family_name, 
                      control=lme4::glmerControl(optimizer = "nloptwrap"), ...)
}

multinomial <- function (link = "logit") {
  linktemp <- substitute(link)
  if (!is.character(linktemp)) 
    linktemp <- deparse(linktemp)
  okLinks <- c("logit", "probit", "cloglog", "cauchit", "log")
  if (linktemp %in% okLinks) 
    stats <- make.link(linktemp)
  else if (is.character(link)) {
    stats <- make.link(link)
    linktemp <- link
  }
  else {
    if (inherits(link, "link-glm")) {
      stats <- link
      if (!is.null(stats$name)) 
        linktemp <- stats$name
    }
    else {
      stop(gettextf("link \"%s\" not available for multinomial family; available links are %s", 
                    linktemp, paste(sQuote(okLinks), collapse = ", ")), 
           domain = NA)
    }
  }
  variance <- function(mu) mu * (1 - mu)
  validmu <- function(mu) all(is.finite(mu)) && all(mu > 0 & mu < 1)
#  dev.resids <- function(y, mu, wt) .Call(stats:::C_binomial_dev_resids, y, mu, wt)
  
  dev.resids <- function(y, mu, wt) .Call(multinomial_dev_resids, y, mu, wt)
  
  aic <- function(y, n, mu, wt, dev) {
    m <- if (any(n > 1)) 
      n
    else wt
    -2 * sum(ifelse(m > 0, (wt/m), 0) * dbinom(round(m * y), round(m), mu, log = TRUE))
  }
 
  initialize <- expression({
    if (NCOL(y) == 1) {
      if (is.factor(y)) y <- y != levels(y)[1L]
      n <- rep.int(1, nobs)
      y[weights == 0] <- 0
      if (any(y < 0)) stop("y values must be 0 <= y <= K")
      #mustart <- (weights * y + 0.5)/(weights + 1)
      # TODO: proper initialization for multinomial
      mustart <- runif(length(y))
      m <- weights * y
      if (any(abs(m - round(m)) > 0.001)) warning("non-integer #successes in a binomial glm!")
    } else if (NCOL(y) == 2) {
      if (any(abs(y - round(y)) > 0.001)) warning("non-integer counts in a binomial glm!")
      n <- y[, 1] + y[, 2]
      y <- ifelse(n == 0, 0, y[, 1]/n)
      weights <- weights * n
      mustart <- (n * y + 0.5)/(n + 1)
    } else stop("for the 'binomial' family, y must be a vector of 0 and 1's\nor a 2 column matrix where col 1 is no. successes and col 2 is no. failures")
  })

  simfun <- function(object, nsim) {
    ftd <- fitted(object)
    n <- length(ftd)
    ntot <- n * nsim
    wts <- object$prior.weights
    if (any(wts%%1 != 0)) 
      stop("cannot simulate from non-integer prior.weights")
    if (!is.null(m <- object$model)) {
      y <- model.response(m)
      if (is.factor(y)) {
        yy <- factor(1 + rbinom(ntot, size = 1, prob = ftd), 
                     labels = levels(y))
        split(yy, rep(seq_len(nsim), each = n))
      }
      else if (is.matrix(y) && ncol(y) == 2) {
        yy <- vector("list", nsim)
        for (i in seq_len(nsim)) {
          Y <- rbinom(n, size = wts, prob = ftd)
          YY <- cbind(Y, wts - Y)
          colnames(YY) <- colnames(y)
          yy[[i]] <- YY
        }
        yy
      }
      else rbinom(ntot, size = wts, prob = ftd)/wts
    }
    else rbinom(ntot, size = wts, prob = ftd)/wts
  }

  linkfun <- stats$linkfun #function(mu) 0
  
  structure(list(family = "multinomial", link = linktemp, linkfun = linkfun, 
                 linkinv = stats$linkinv, variance = variance, dev.resids = dev.resids, 
                 aic = aic, mu.eta = stats$mu.eta, initialize = initialize, 
                 validmu = validmu, valideta = stats$valideta, simulate = simfun), 
            class = "family")
  
  # structure(list(family = "multinomial", link = linktemp, linkfun = stats$linkfun, 
  #                linkinv = stats$linkinv, variance = variance, dev.resids = dev.resids, 
  #                aic = aic, mu.eta = stats$mu.eta, initialize = initialize, 
  #                validmu = validmu, valideta = stats$valideta, simulate = simfun), 
  #           class = "family")
}
