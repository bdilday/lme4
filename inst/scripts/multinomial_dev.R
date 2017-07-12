
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
  df1$outcome <- mm
  df1
}

make_multinomial_df <- function(df1) {
  column_names <- df1 %>% 
    subset(select = !grepl('outcome', names(df1))) %>% names()
  all_outcomes <- unique(df1$outcome)
  tmp <- data_frame(outcome=df1$outcome)
  nl <- length(all_outcomes) - 1
  for (k in 1:nl) {
    for (column_name in column_names) {
      new_column_name <- paste(column_name, k, sep='.')
      tmp[[new_column_name]] <- as.factor(df1[[column_name]])
    }
  }
  
  tmp
}

#' @export 
step_through_binomial <- function(nlim=1000) {

}

#' @export
test_multinomial_data <- function(nlim=NULL) {
  BProDRA::generate_model_df(nlim=nlim, year=2016)  
}

#' @export
test_multinomial_data_X <- function(nlim=NULL) {
  column_names <- c("BAT_ID", "PIT_ID", "HOME_TEAM_ID", "bid", "pid", "sid")
  model_df <- BProDRA::generate_model_df(nlim=nlim, year=2016)
  all_outcomes <- unique(model_df$ev$outcome) %>% sort()
  ref_outcome <- all_outcomes[[1]] 
  nl <- length(all_outcomes)
  ev <- list()
  
  for (k in 2:nl) {
    cc <- which(model_df$ev$outcome == all_outcomes[[k]])
    tmp <- model_df$ev[cc,]
    for (column_name in column_names) {
      tmp[[column_name]] <- paste(tmp[[column_name]], k, sep='_')
      tmp$outcome <- 1
    }
    ev[[k]] <- tmp
    
    cc <- which(model_df$ev$outcome == all_outcomes[[1]])
    tmp <- model_df$ev[cc,]
    for (column_name in column_names) {
      tmp[[column_name]] <- paste(tmp[[column_name]], k, sep='_')
      tmp$outcome <- 0
    }
    ev[[100 + k]] <- tmp
  }

  tmp <- purrr::reduce(ev, rbind.data.frame)
  model_df$ev <- tmp
  rm(tmp)
  model_df
}

do_multinomial_as_retrms <- function(nlim=NULL) {
#  column_names <- c("BAT_ID", "PIT_ID", "HOME_TEAM_ID")
  column_names <- c("BAT_ID")
  model_df <- BProDRA::generate_model_df(nlim=nlim, year=2016)
  all_outcomes <- unique(model_df$ev$outcome) %>% sort()
  nl <- length(all_outcomes)
  ev <- list()
  
  tmp <- model_df$ev  
  for (k in 1:nl) {
    for (column_name in column_names) {
      new_column_name <- paste(column_name, k, sep='.')
      tmp[[new_column_name]] <- as.factor(tmp[[column_name]])
    }
  }
  
  tmp %<>% mutate(idcase=row_number())  
  glf <- lme4::glFormula(
    outcome ~ (1|BAT_ID.1) + 
      (1|BAT_ID.2) + 
      (1|BAT_ID.3) + 
      (1|BAT_ID.4), 
    data=tmp, family=binomial, nAGQ=0
    )
 
   ref_index <- last(all_outcomes)
   cc = which(tmp$outcome == ref_index)
   tmp[cc,]$outcome <- 0
   
   for (idx in 1:(ref_index-1)) {
     cc <- which(tmp$outcome == idx)
     
   }
   
}

#' @export
test_multinomial_data_mlogit <- function(nlim=NULL) {
  column_names <- c("BAT_ID", "PIT_ID", "HOME_TEAM_ID")
  model_df <- BProDRA::generate_model_df(nlim=nlim, year=2016)
  all_outcomes <- unique(model_df$ev$outcome) %>% sort()
  nl <- length(all_outcomes)
  ev <- list()

  tmp <- model_df$ev  
  for (k in 1:nl) {
    for (column_name in column_names) {
      new_column_name <- paste(column_name, k, sep='.')
      tmp[[new_column_name]] <- as.factor(tmp[[column_name]])
    }
  }

  tmp %<>% mutate(idcase=row_number())
  model_df$ev <- tmp
  rm(tmp)
  model_df
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
      mustart <- y/NCOL(y)
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

  # TODO: finish this
  linkfun <- function(mu) {
    dd <- dim(mu)
    nr <- dim[1]
    z <- mu
  }
  
  linkinv <- function() {
    1
  }
  
  mu.eta <- function() {
    1
  } 
  
  valideta <- function() {
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
  devfun <- do.call(mkGlmerDevfun, c(glmod, list(verbose = verbose,
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
