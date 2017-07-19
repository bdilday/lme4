
#' simulated multinomial data
#' 
#' simulates data with 3 x 2 factors
#' 
#' @export
simulate_data <- function(rseed=101, n_matchup=10) {
  set.seed(rseed)
  sim_coefs <- list('A' = c(1, 1, 0), 
                    'B' = c(-1, 1, 0), 
                    'C' = c(0, 0, 1)
  )
  
  sim_entry <- function(f1, f2) {
    eta <- sim_coefs[[f1]] + sim_coefs[[f2]]
    multinom_probs <- exp(eta)/sum(exp(eta))
    multinom_draw <- as.vector(rmultinom(1, 1, multinom_probs))
    which(multinom_draw > 0)
  }
  list_of_factors <- names(sim_coefs)
  factor1 <- rep(list_of_factors, each=n_matchup * length(list_of_factors))
  factor2 <- rep(lapply(list_of_factors, rep, n_matchup), length(list_of_factors)) %>% unlist()

  df1 <- data_frame(fB=factor1, fP=factor2)
  results_vector <- sapply(1:nrow(df1), function(idx) {
    r <- df1[idx,]
    sim_entry(r$fB, r$fP)
  })
  df1$outcome <- results_vector
  df1
}

#' simulated multinomial data
#' 
#' simulates data with 3 x 2 factors
#' 
#' @export
simulate_data_as_matrix_resp <- function(rseed=101, n_matchup=10) {
  set.seed(rseed)
  sim_coefs <- list('A' = c( 1, 1, 0), 
                    'B' = c(-1, 1, 0), 
                    'C' = c( 0, 0, 1)
  )
  
  sim_entry <- function(f1, f2) {
    eta <- sim_coefs[[f1]] + sim_coefs[[f2]]
    multinom_probs <- exp(eta)/sum(exp(eta))
    multinom_draw <- as.vector(rmultinom(1, 1, multinom_probs))
    which(multinom_draw > 0)
  }
  list_of_factors <- names(sim_coefs)
  factor1 <- rep(list_of_factors, each=n_matchup * length(list_of_factors))
  factor2 <- rep(lapply(list_of_factors, rep, n_matchup), length(list_of_factors)) %>% unlist()
  
  df1 <- data_frame(fB=factor1, fP=factor2)
  results_vector <- sapply(1:nrow(df1), function(idx) {
    r <- df1[idx,]
    sim_entry(r$fB, r$fP)
  })
  K_class <- length(sim_coefs[[1]])
  mm <- matrix(rep(0, (K_class) * nrow(df1)), ncol=(K_class))
  for(i in 1:nrow(df1)) {mm[i,results_vector[i]] <- 1}
  df1$outcome <- mm[,1:(K_class-1)]
  df1
}

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


#' @export
multinomial <- function (link = "multiclasslogit") {
  linktemp <- substitute(link)
  if (!is.character(linktemp)) 
    linktemp <- deparse(linktemp)
  okLinks <- c("multiclasslogit", "multiclassprobit")
  stopifnot(linktemp %in% okLinks) 
  
  variance <- function(mu) {
    mu * (1 - mu)
  }
  validmu <- function(mu) {
    all(is.finite(mu)) && all(mu > 0 & mu < 1)
  }
#  dev.resids <- function(y, mu, wt) .Call(stats:::C_binomial_dev_resids, y, mu, wt)
  
  dev.resids <- function(y, mu, wt) {
    .Call(multinomial_dev_resids, y, mu, wt)
    #1
   # .Call(stats:::C_binomial_dev_resids, y, mu, wt)
  }
  
  aic <- function(y, n, mu, wt, dev) {
    m <- if (any(n > 1)) 
      n
    else wt
    -2 * sum(ifelse(m > 0, (wt/m), 0) * dbinom(round(m * y), round(m), mu, log = TRUE))
  }
 
  initialize <- expression({
    if (NCOL(y) == 1) {
      stop("for the 'multinomial' family, y must be an N x K matrix of 0's and 1's")
    } else { 
      if (any(abs(y - round(y)) > 0.001)) warning("non-integer counts in a multinomial glm!")
      mustart <- matrix(1/(1+NCOL(y)), nrow=dim(y)[1], ncol=dim(y)[2])
    }
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

  linkfun <- function(mu) {
    dd <- dim(mu)
    stopifnot(!is.null(dd))
    stopifnot(length(dd) == 2)
    z <- mu/(1 - rowSums(mu))
    log(z)
  }

  # TODO: finish this  
  linkinv <- function(eta) {
    dd <- dim(eta)
    message('dim(eta) ', dd, ' ', is.null(dd), ' ', length(eta))
    stopifnot(!is.null(dd))
    stopifnot(length(dd) == 2)
    exp(eta)/(1+rowSums(exp(eta)))
  }
  
  mu.eta <- function(eta) {
    mu <- linkinv(eta)
    mu * (1 - mu)
  } 
  
  valideta <- function(eta) {
    TRUE
  }
  
  stats <- structure(
    list(linkfun = linkfun, 
         linkinv = linkinv, 
         mu.eta = mu.eta,
         valideta = valideta, 
         name = linktemp), 
    class = "link-glm")
  
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


#' @export
##' Fit a generalized linear mixed model (GLMM)
glmer_vec <- function(formula, data=NULL, family = gaussian,
                  control = glmerControl(), start = NULL, verbose = 0L, nAGQ = 1L,
                  subset, weights, na.action, offset,
                  contrasts = NULL, mustart, etastart, devFunOnly = FALSE, ...)
{
  if (!inherits(control, "glmerControl")) {
    if(!is.list(control)) stop("'control' is not a list; use glmerControl()")
    ## back-compatibility kluge
    if (class(control)[1]=="lmerControl") {
      warning("please use glmerControl() instead of lmerControl()",
              immediate.=TRUE)
      control <-
        ## unpack sub-lists
        c(control[!names(control) %in% c("checkConv","checkControl")],
          control$checkControl,control$checkConv)
      control["restart_edge"] <- NULL ## not implemented for glmer
    } else {
      msg <- "Use control=glmerControl(..) instead of passing a list"
      if(length(cl <- class(control))) {
        msg <- paste(msg, "of class", dQuote(cl[1]))
      }
      warning(msg, immediate.=TRUE)
    }
    
    control <- do.call(glmerControl, control)
  }
  mc <- mcout <- match.call()
  
  ## family-checking code duplicated here and in glFormula (for now) since
  ## we really need to redirect at this point; eventually deprecate formally
  ## and clean up
  if (is.character(family))
    family <- get(family, mode = "function", envir = parent.frame(2))
  if( is.function(family)) family <- family()
  if (isTRUE(all.equal(family, gaussian()))) {
    ## redirect to lmer (with warning)
    warning("calling glmer() with family=gaussian (identity link) as a shortcut to lmer() is deprecated;",
            " please call lmer() directly")
    mc[[1]] <- quote(lme4::lmer)
    mc["family"] <- NULL            # to avoid an infinite loop
    return(eval(mc, parent.frame()))
  }
  
  ## see https://github.com/lme4/lme4/issues/50
  ## parse the formula and data
  mc[[1]] <- quote(lme4::glFormula)
  glmod <- eval(mc, parent.frame(1L))
  mcout$formula <- glmod$formula
  glmod$formula <- NULL
  
  ## create deviance function for covariance parameters (theta)
  
  nAGQinit <- if(control$nAGQ0initStep) 0L else 1L
  devfun <- do.call(mkGlmerVecDevfun, c(glmod, list(verbose = verbose,
                                                    control = control,
                                                    nAGQ = nAGQinit)))
  if (nAGQ==0 && devFunOnly) return(devfun)
  ## optimize deviance function over covariance parameters
  
  ## FIXME: perhaps should be in glFormula instead??
  if (is.list(start)) {
    start.bad <- setdiff(names(start),c("theta","fixef"))
    if (length(start.bad)>0) {
      stop(sprintf("bad name(s) for start vector (%s); should be %s and/or %s",
                   paste(start.bad,collapse=", "),
                   shQuote("theta"),
                   shQuote("fixef")),call.=FALSE)
    }
    if (!is.null(start$fixef) && nAGQ==0)
      stop("should not specify both start$fixef and nAGQ==0")
  }
  
  ## FIX ME: allow calc.derivs, use.last.params etc. if nAGQ=0
  if(control$nAGQ0initStep) {
    opt <- optimizeGlmer(devfun,
                         optimizer = control$optimizer[[1]],
                         ## DON'T try fancy edge tricks unless nAGQ=0 explicitly set
                         restart_edge=if (nAGQ==0) control$restart_edge else FALSE,
                         boundary.tol=if (nAGQ==0) control$boundary.tol else 0,
                         control = control$optCtrl,
                         start=start,
                         nAGQ = 0,
                         verbose=verbose,
                         calc.derivs=FALSE)
  }
  
  if(nAGQ > 0L) {
    
    
    ## update deviance function to include fixed effects as inputs
    devfun <- updateGlmerDevfun(devfun, glmod$reTrms, nAGQ = nAGQ)
    
    if (control$nAGQ0initStep) {
      start <- updateStart(start,theta=opt$par)
    }
    ## if nAGQ0 was skipped
    ## we don't actually need to do anything here, it seems --
    ## getStart gets called again in optimizeGlmer ...
    
    if (devFunOnly) return(devfun)
    ## reoptimize deviance function over covariance parameters and fixed effects
    opt <- optimizeGlmer(devfun,
                         optimizer = control$optimizer[[2]],
                         restart_edge=control$restart_edge,
                         boundary.tol=control$boundary.tol,
                         control = control$optCtrl,
                         start=start,
                         nAGQ=nAGQ,
                         verbose = verbose,
                         stage=2,
                         calc.derivs=control$calc.derivs,
                         use.last.params=control$use.last.params)
  }
  cc <- if (!control$calc.derivs) NULL else {
    if (verbose > 10) cat("checking convergence\n")
    checkConv(attr(opt,"derivs"),opt$par,
              ctrl = control$checkConv,
              lbound=environment(devfun)$lower)
  }
  
  ## prepare output
  mkMerMod(environment(devfun), opt, glmod$reTrms, fr = glmod$fr,
           mc = mcout, lme4conv=cc)
  
}## {glmer}




##' Create an lmerResp, glmResp or nlsResp instance
##'
##' @title Create an lmerResp, glmResp or nlsResp instance
##' @param fr a model frame
##' @param REML logical scalar, value of REML for an lmerResp instance
##' @param family the optional glm family (glmResp only)
##' @param nlenv the nonlinear model evaluation environment (nlsResp only)
##' @param nlmod the nonlinear model function (nlsResp only)
##' @param ... where to look for response information if \code{fr} is missing.
##'   Can contain a model response, \code{y}, offset, \code{offset}, and weights,
##'   \code{weights}.
##' @return an lmerResp or glmResp or nlsResp instance
##' @family utilities
##' @export
mkRespVecMod <- function(fr, REML=NULL, family = NULL, nlenv = NULL, nlmod = NULL, ...)
{
  if(!missing(fr)) {
    y <- model.response(fr)
    offset <- model.offset(fr)
    weights <- model.weights(fr)
    N <- n <- nrow(fr)
    etastart_update <- model.extract(fr, "etastart")
    mustart_update <- model.extract(fr, "mustart")
  } else {
    fr <- list(...)
    y <- fr$y
    N <- n <- NROW(y)
    offset <- fr$offset
    weights <- fr$weights
    etastart_update <- fr$etastart
    mustart_update <- fr$mustart
  }
  if(length(dim(y)) == 1L)
    y <- drop(y) ## avoid problems with 1D arrays and keep names
  
  if(isGLMM <- !is.null(family))
    stopifnot(inherits(family, "family"))
  ## FIXME: may need to add X, or pass it somehow, if we want to use glm.fit
  
  ## test for non-numeric response here to avoid later
  ## confusing error messages from deeper machinery
  if (!is.null(y)) { ## 'y' may be NULL if we're doing simulation
    if(!(is.numeric(y) ||
         ((is.binom <- isGLMM && family$family == "binomial") &&
          (is.factor(y) || is.logical(y))))) {
      if (is.binom)
        stop("response must be numeric or factor")
      else {
        if (is.logical(y))
          y <- as.integer(y)
        else stop("response must be numeric")
      }
    }
    if(!all(is.finite(y)))
      stop("NA/NaN/Inf in 'y'") # same msg as from lm.fit()
  }
  
  rho <- new.env()
  rho$y_vec <- if (is.null(y)) numeric(0) else y
  if (!is.null(REML)) rho$REML <- REML
  rho$etastart <- etastart_update
  rho$mustart <- mustart_update
  rho$start <- attr(fr,"start")
  if (!is.null(nlenv)) {
    stopifnot(is.language(nlmod),
              is.environment(nlenv),
              is.numeric(val <- eval(nlmod, nlenv)),
              length(val) == n,
              ## FIXME?  Restriction, not present in ole' nlme():
              is.matrix(gr <- attr(val, "gradient")),
              is.numeric(gr),
              nrow(gr) == n,
              !is.null(pnames <- colnames(gr)))
    N <- length(gr)
    rho$mu <- as.vector(val)
    rho$sqrtXwt <- as.vector(gr)
    rho$gam <- ## FIXME more efficient  mget(pnames, envir=nlenv)
      unname(unlist(lapply(pnames,
                           function(nm) get(nm, envir=nlenv))))
  }
  rho$offset <- if (!is.null(offset)) {
    if (length(offset) == 1L) offset <- rep.int(offset, N)
    else stopifnot(length(offset) == N)
    unname(offset)
  } else rep.int(0, N)
  rho$weights <- if (!is.null(weights)) {
    stopifnot(length(weights) == n, all(weights >= 0))
    unname(weights)
  } else rep.int(1, n)
  
  if(isGLMM) {
    ## need weights for initializing evaluation
    rho$nobs <- n
    ## allow trivial objects, e.g. for simulation
    if (length(y)>0) eval(family$initialize, rho)
    ## ugh. this *is* necessary;
    ##  family$initialize *ignores* mustart in env, overwrites!
    ## see ll 180-182 of src/library/stats/R/glm.R
    ## https://github.com/wch/r-source/search?utf8=%E2%9C%93&q=mukeep
    if (!is.null(mustart_update)) rho$mustart <- mustart_update
    ## family$initialize <- NULL     # remove clutter from str output
    ll <- as.list(rho)
    rho$y_vec <- y
    rho$mu_vec <- rho$mustart
    rho$eta_vec <- rho$etastart
    
    ll[['initial_dim']] <- dim(ll$y_vec)
    ll[['k_class']] <- dim(ll$y_vec)[[2]]
    ll[['calling_source']] <- 'mkRespVecMod'
    ans <- do.call(new, c(list(Class="glmVecResp", family=family),
                          ll[setdiff(names(ll), c("m", "nobs", "mustart"))]))
    if (length(y)>0) {
      es <- etastart_update
      if (!is.null(es)) {
        upmu_arg <- es 
      } else {
        upmu_arg <- family$linkfun(rho$mustart)
      }
      ans$updateMu(upmu_arg)
    }
    
    ans
    
  } else if (is.null(nlenv)) ## lmer
    do.call(lmerResp$new, as.list(rho))
  else ## nlmer
    do.call(nlsResp$new,
            c(list(nlenv=nlenv,
                   nlmod=substitute(~foo, list(foo=nlmod)),
                   pnames=pnames), as.list(rho)))
}

#' @export
fix_multinomial_environment <- function(rho) {
  nn <- names(rho)
  x <- 1
}

#' @export 
step_through_multinomial <- function(df1=NULL, n_matchup=10, verbose=0) {
  if (is.null(df1)) {
    df1 <- simulate_data_as_matrix_resp(n_matchup = 10)    
  }
  
  mc <- mcout <- match.call()
  control <- lme4::glmerControl(optimizer = "nloptwrap")
  nAGQinit <- 0
  
  glf <- lme4::glVecFormula(outcome ~ (1|fB) + (1|fP), 
                            data=df1, nAGQ = 0, 
                            family = multinomial, 
                            control=control, verbose=1000)
  
  mcout$formula <- glf$formula
  glf$formula <-  NULL # ??
  
  devfun <- do.call(mkGlmerVecDevfun, c(glf, list(verbose = verbose,
                                                  control = control,
                                                  nAGQ = nAGQinit)))
  
  opt <- optimizeGlmer(devfun,
                       optimizer = control$optimizer[[1]],
                       ## DON'T try fancy edge tricks unless nAGQ=0 explicitly set
                       restart_edge=if (nAGQ==0) control$restart_edge else FALSE,
                       boundary.tol=if (nAGQ==0) control$boundary.tol else 0,
                       control = control$optCtrl,
                       start=start,
                       nAGQ = 0,
                       verbose=verbose,
                       calc.derivs=FALSE)
  
  out <- mkMerMod(environment(devfun), opt, glmod$reTrms, fr = glmod$fr,
           mc = mcout, lme4conv=cc)
  
  out
}


