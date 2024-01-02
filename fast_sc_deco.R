

library(GJRM)



sigmoid <- function(x){
  exp(x)/(exp(x)+1)
}

logit <- function(x){
  log(x/(1-x))
}





gjrm <- function(formula, data = list(), weights = NULL, subset = NULL, 
                 copula = "N", copula2 = "N", margins, model, dof = 3, dof2 = 3, 
                 ordinal = FALSE, surv = FALSE, cens1 = NULL, cens2 = NULL, 
                 cens3 = NULL, dep.cens = FALSE, upperBt1 = NULL, upperBt2 = NULL, 
                 gamlssfit = FALSE, fp = FALSE, infl.fac = 1, rinit = 1, 
                 rmax = 1000,  iterlimsp = 50, tolsp = 1e-07, gc.l = FALSE, 
                 parscale, extra.regI = "t", k1.tvc = 0, k2.tvc = 0, knots = NULL, 
                 penCor = "unpen", sp.penCor = 3, Chol = FALSE, gamma = 1, 
                 w.alasso = NULL, drop.unused.levels = TRUE, min.dn = 1e-40, 
                 min.pr = 1e-16, max.pr = 0.999999) 
{
  BivD <- copula
  BivD2 <- copula2
  Model <- model
  if (missing(margins)) 
    stop("You must choose the margins' values.")
  if (missing(Model)) 
    stop("You must choose a model type.")
  if (margins[1] == "PH" && surv == TRUE) 
    margins[1] <- "cloglog"
  if (margins[1] == "PO" && surv == TRUE) 
    margins[1] <- "logit"
  if (margins[2] == "PH" && surv == TRUE) 
    margins[2] <- "cloglog"
  if (margins[2] == "PO" && surv == TRUE) 
    margins[2] <- "logit"
  bl <- c("probit", "logit", "cloglog")
  v.rB1 <- upperBt1
  v.rB2 <- upperBt2
  if (Model == "ROY") {
    L <- eval(substitute(SemiParROY(formula, data, weights, 
                                    subset, BivD1 = BivD, BivD2, margins, dof1 = dof, 
                                    dof2, gamlssfit, fp, infl.fac, rinit, rmax, iterlimsp, 
                                    tolsp, gc.l, parscale, extra.regI, knots = knots, 
                                    drop.unused.levels = drop.unused.levels, min.dn = min.dn, 
                                    min.pr = min.pr, max.pr = max.pr), list(weights = weights)))
  }
  else {
    if (surv == FALSE && ordinal == FALSE) {
      if ((margins[1] %in% bl && margins[2] %in% bl && 
           is.na(margins[3])) || (margins[1] %in% bl && 
                                  !(margins[2] %in% bl) && Model == "B" && is.na(margins[3]))) {
        L <- eval(substitute(SemiParBIV(formula, data, 
                                        weights, subset, Model, BivD, margins, dof, 
                                        gamlssfit, fp, hess = TRUE, infl.fac, rinit, 
                                        rmax, iterlimsp, tolsp, gc.l, parscale, extra.regI, 
                                        intf = TRUE, theta.fx = NULL, knots = knots, 
                                        drop.unused.levels = drop.unused.levels, min.dn = min.dn, 
                                        min.pr = min.pr, max.pr = max.pr), list(weights = weights)))
      }
    }
    if (surv == FALSE && ordinal == TRUE) {
      if ((margins[1] %in% bl && margins[2] %in% bl && 
           is.na(margins[3])) || (margins[1] %in% bl && 
                                  !(margins[2] %in% bl) && is.na(margins[3]))) {
        L <- eval(substitute(CopulaCLM(formula, data, 
                                       weights, subset, Model, BivD, margins, dof, 
                                       gamlssfit, fp, hess = TRUE, infl.fac, rinit, 
                                       rmax, iterlimsp, tolsp, gc.l, parscale, extra.regI, 
                                       intf = TRUE, theta.fx = NULL, knots = knots, 
                                       drop.unused.levels = drop.unused.levels, min.dn = min.dn, 
                                       min.pr = min.pr, max.pr = max.pr), list(weights = weights)))
      }
    }
    if (margins[1] %in% bl && !(margins[2] %in% bl) && surv == 
        FALSE && is.na(margins[3]) && Model == "BSS" && 
        ordinal == FALSE) {
      L <- eval(substitute(copulaSampleSel(formula, data, 
                                           weights, subset, BivD, margins, dof, fp, infl.fac, 
                                           rinit, rmax, iterlimsp, tolsp, gc.l, parscale, 
                                           extra.regI, knots, drop.unused.levels = drop.unused.levels, 
                                           min.dn = min.dn, min.pr = min.pr, max.pr = max.pr), 
                           list(weights = weights)))
    }
    if (!is.na(margins[3])) {
      if (margins[1] %in% bl && margins[2] %in% bl && 
          margins[3] %in% bl && surv == FALSE && ordinal == 
          FALSE) {
        L <- eval(substitute(SemiParTRIV(formula, data, 
                                         weights, subset, Model, margins, penCor, sp.penCor, 
                                         approx = FALSE, Chol, infl.fac, gamma, w.alasso, 
                                         rinit, rmax, iterlimsp, tolsp, gc.l, parscale, 
                                         extra.regI, knots, drop.unused.levels = drop.unused.levels, 
                                         min.dn = min.dn, min.pr = min.pr, max.pr = max.pr), 
                             list(weights = weights)))
      }
    }
    if ((!(margins[1] %in% bl) || surv == TRUE) && ordinal == 
        FALSE) {
      robust <- FALSE
      t.c = 3
      sp <- qu.mag <- y1.y2 <- y1.cy2 <- cy1.y2 <- cy1.cy2 <- cy <- cy1 <- gamlss1 <- gamlss2 <- gam1 <- gam2 <- y1m <- y2m <- indexTeq1 <- indexTeq2 <- NULL
      i.rho <- log.sig2.2 <- log.nu.2 <- log.nu.1 <- log.sig2.1 <- dof.st <- NULL
      end <- X3.d2 <- X4.d2 <- X5.d2 <- X6.d2 <- X7.d2 <- X8.d2 <- l.sp3 <- l.sp4 <- l.sp5 <- l.sp6 <- l.sp7 <- l.sp8 <- l.sp9 <- 0
      sp1 <- sp2 <- NULL
      sp3 <- gp3 <- gam3 <- X3 <- sp4 <- gp4 <- gam4 <- X4 <- sp5 <- gp5 <- gam5 <- X5 <- gam9 <- NULL
      sp6 <- gp6 <- gam6 <- X6 <- sp7 <- gp7 <- gam7 <- X7 <- sp8 <- gp8 <- gam8 <- X8 <- sp9 <- NULL
      c11 <- c10 <- c01 <- c00 <- NA
      cens1Mix <- cens2Mix <- NULL
      Sl.sf <- NULL
      sp.method <- "perf"
      Xd1 <- Xd2 <- mono.sm.pos1 <- mono.sm.pos2 <- mono.sm.pos <- NULL
      surv.flex <- FALSE
      Deq1 <- pos.pbeq1 <- Deq2 <- pos.pbeq2 <- list()
      BivD2 <- c("C0C90", "C0C270", "C180C90", "C180C270", 
                 "J0J90", "J0J270", "J180J90", "J180J270", "G0G90", 
                 "G0G270", "G180G90", "G180G270", "GAL0GAL90", 
                 "GAL0GAL270", "GAL180GAL90", "GAL180GAL270")
      opc <- c("N", "C0", "C90", "C180", "C270", "J0", 
               "J90", "J180", "J270", "G0", "G90", "G180", 
               "G270", "F", "AMH", "FGM", "T", "PL", "HO", 
               "GAL0", "GAL90", "GAL180", "GAL270")
      scc <- c("C0", "C180", "GAL0", "GAL180", "J0", "J180", 
               "G0", "G180", BivD2)
      sccn <- c("C90", "C270", "GAL90", "GAL270", "J90", 
                "J270", "G90", "G270")
      m2 <- c("N", "GU", "rGU", "LO", "LN", "WEI", "iG", 
              "GA", "BE", "FISK", "GP", "GPII", "GPo")
      m3 <- c("DAGUM", "SM", "TW")
      m1d <- c("PO", "ZTP", "DGP0")
      m2d <- c("NBI", "NBII", "PIG", "DGP", "DGPII")
      m3d <- c("DEL", "SICHEL")
      ct <- data.frame(c(opc), c(1:14, 55, 56, 57, 60, 
                                 61, 62:65))
      cta <- data.frame(c(opc), c(1, 3, 23, 13, 33, 6, 
                                  26, 16, 36, 4, 24, 14, 34, 5, 55, 56, 2, 60, 
                                  61, 62:65))
      if (BivD %in% BivD2) {
        if (BivD %in% BivD2[1:4]) 
          BivDt <- "C0"
        if (BivD %in% BivD2[5:12]) 
          BivDt <- "J0"
        if (BivD %in% BivD2[13:16]) 
          BivDt <- "C0"
        nC <- ct[which(ct[, 1] == BivDt), 2]
        nCa <- cta[which(cta[, 1] == BivDt), 2]
      }
      if (!(BivD %in% BivD2)) {
        nC <- ct[which(ct[, 1] == BivD), 2]
        nCa <- cta[which(cta[, 1] == BivD), 2]
      }
      if (!is.list(formula)) 
        stop("You must specify a list of equations.")
      l.flist <- length(formula)
      form.check(formula, l.flist)
      cl <- match.call()
      mf <- match.call(expand.dots = FALSE)
      pred.varR <- pred.var(formula, l.flist)
      v1 <- pred.varR$v1
      v2 <- pred.varR$v2
      pred.n <- pred.varR$pred.n
      if (!is.null(v.rB1)) 
        pred.n <- c(pred.n, v.rB1)
      if (!is.null(v.rB2)) 
        pred.n <- c(pred.n, v.rB2)
      fake.formula <- paste(v1[1], "~", paste(pred.n, 
                                              collapse = " + "))
      environment(fake.formula) <- environment(formula[[1]])
      mf$formula <- fake.formula
      mf$upperBt1 <- mf$upperBt2 <- mf$min.dn <- mf$min.pr <- mf$max.pr <- mf$dep.cens <- mf$ordinal <- mf$Model <- mf$model <- mf$knots <- mf$k1.tvc <- mf$k2.tvc <- mf$surv <- mf$BivD <- mf$copula <- mf$copula2 <- mf$margins <- mf$fp <- mf$dof <- mf$infl.fac <- mf$rinit <- mf$rmax <- mf$iterlimsp <- mf$tolsp <- mf$gc.l <- mf$parscale <- mf$extra.regI <- mf$gamlssfit <- NULL
      mf$drop.unused.levels <- drop.unused.levels
      if (surv == TRUE) 
        mf$na.action <- na.pass
      mf[[1]] <- as.name("model.frame")
      data <- eval(mf, parent.frame())
      if (surv == TRUE) {
        if (!("(cens1)" %in% names(data)) && margins[1] %in% 
            bl) 
          stop("You must provide both censoring indicators.")
        if (!("(cens2)" %in% names(data)) && margins[2] %in% 
            bl) 
          stop("You must provide both censoring indicators.")
      }
      if (gc.l == TRUE) 
        gc()
      if (!("(weights)" %in% names(data))) {
        weights <- rep(1, dim(data)[1])
        data$weights <- weights
        names(data)[length(names(data))] <- "(weights)"
      }
      else weights <- data[, "(weights)"]
      if (!("(cens1)" %in% names(data))) {
        cens1 <- rep(0, dim(data)[1])
        data$cens1 <- cens1
        names(data)[length(names(data))] <- "(cens1)"
      }
      else cens1 <- data[, "(cens1)"]
      if (!("(cens2)" %in% names(data))) {
        cens2 <- rep(0, dim(data)[1])
        data$cens2 <- cens2
        names(data)[length(names(data))] <- "(cens2)"
      }
      else cens2 <- data[, "(cens2)"]
      if (!("(cens3)" %in% names(data))) {
        cens3 <- rep(0, dim(data)[1])
        data$cens3 <- cens3
        names(data)[length(names(data))] <- "(cens3)"
      }
      else cens3 <- data[, "(cens3)"]
      if (surv == TRUE) {
        if (is.factor(cens1) && !is.factor(cens2)) 
          stop("Both censoring indicators have to be factor variables for mixed censoring case.")
        if (!is.factor(cens1) && is.factor(cens2)) 
          stop("Both censoring indicators have to be factor variables for mixed censoring case.")
      }
      if (surv == TRUE && is.factor(cens1) && !is.null(v.rB1)) 
        data[!(cens1 == "I"), v.rB1] <- data[!(cens1 == 
                                                 "I"), v1[1]]
      if (surv == TRUE && is.factor(cens2) && !is.null(v.rB2)) 
        data[!(cens2 == "I"), v.rB2] <- data[!(cens2 == 
                                                 "I"), v2[1]]
      if (surv == TRUE) {
        if (any(is.na(data[, v1[1]]) | is.na(data[, 
                                                  v2[1]]))) 
          stop("Time to event with NA's. Please check your time covariates.")
        actual.NAs = as.numeric(which(apply(apply(data, 
                                                  1, is.na), 2, any)))
        data <- na.omit(data)
        if (length(actual.NAs) > 0) {
          cens1 <- cens1[-actual.NAs]
          cens2 <- cens2[-actual.NAs]
        }
      }
      n <- dim(data)[1]
      if (surv == TRUE && is.factor(cens1) && is.factor(cens2)) {
        cens1Mix <- cens1
        cens2Mix <- cens2
        cens1 <- cens2 <- rep(1, n)
      }
      M <- list(m1d = m1d, m2 = m2, m2d = m2d, m3 = m3, 
                m3d = m3d, BivD = BivD, bl = bl, robust = robust, 
                opc = opc, extra.regI = extra.regI, margins = margins, 
                BivD2 = BivD2, dof = dof, surv = surv, c1 = cens1, 
                c2 = cens2, c3 = cens3, dep.cens = dep.cens)
      M$K1 <- NULL
      M$type.cens1 <- M$type.cens2 <- "R"
      pream.wm(formula, margins, M, l.flist)
      formula.eq1 <- formula[[1]]
      formula.eq2 <- formula[[2]]
      form.eq12R <- form.eq12(formula.eq1, data, v1, margins[1], 
                              m1d, m2d)
      formula.eq1 <- form.eq12R$formula.eq1
      formula.eq1r <- form.eq12R$formula.eq1r
      y1 <- form.eq12R$y1
      y1.test <- form.eq12R$y1.test
      y1m <- form.eq12R$y1m
      if (surv == FALSE) 
        gam1 <- eval(substitute(gam(formula.eq1, gamma = infl.fac, 
                                    weights = weights, data = data, knots = knots, 
                                    drop.unused.levels = drop.unused.levels), 
                                list(weights = weights)))
      if (surv == TRUE && margins[1] %in% c(m2, m3) && 
          margins[2] %in% bl) 
        gam1 <- eval(substitute(gam(formula.eq1, gamma = infl.fac, 
                                    weights = weights, data = data, knots = knots, 
                                    drop.unused.levels = drop.unused.levels), 
                                list(weights = weights)))
      else {
        if (surv == TRUE && !(margins[1] %in% bl)) 
          gam1 <- eval(substitute(gam(formula.eq1, gamma = infl.fac, 
                                      weights = weights * cens1, data = data, 
                                      knots = knots, drop.unused.levels = drop.unused.levels), 
                                  list(weights = weights, cens1 = cens1)))
      }
      if (surv == TRUE && margins[1] %in% bl) {
        surv.flex <- TRUE
        f.eq1 <- form.eq12R$f.eq1
        data$urcfcphmwicu <- seq(-10, 10, length.out = dim(data)[1])
        tempb <- eval(substitute(gam(f.eq1, family = cox.ph(), 
                                     data = data, weights = cens1, drop.unused.levels = drop.unused.levels), 
                                 list(cens1 = cens1)))
        data$Sh <- as.vector(mm(predict(tempb, type = "response"), 
                                min.pr = min.pr, max.pr = max.pr))
        cens11 <- ifelse(cens1 == 0, 1e-07, cens1)
        gam1 <- eval(substitute(scam(formula.eq1, gamma = infl.fac, 
                                     weights = weights * cens11, data = data), 
                                list(weights = weights, cens11 = cens11)))
        lsgam1 <- length(gam1$smooth)
        if (lsgam1 == 0) 
          stop("You must use at least a monotonic smooth function of time in the first equation.")
        clsm <- ggr <- NA
        for (i in 1:lsgam1) {
          clsm[i] <- class(gam1$smooth[[i]])[1]
        }
        if (sum(as.numeric(clsm[1] %in% c("mpi.smooth"))) == 
            0) 
          stop("You must have a monotonic smooth of time and it has to be the first to be included.")
        l.sp1 <- length(gam1$sp)
        if (l.sp1 != 0) 
          sp1 <- gam1$sp
        sp1[1] <- 1
        gam.call <- gam1$call
        gam.call$sp <- sp1
        gam1 <- eval(gam.call)
        j <- 1
        for (i in 1:lsgam1) {
          if (max(as.numeric(grepl(v1[1], gam1$smooth[[i]]$term))) != 
              0 && clsm[i] == "mpi.smooth") 
            mono.sm.pos1 <- c(mono.sm.pos1, c(gam1$smooth[[i]]$first.para:gam1$smooth[[i]]$last.para))
        }
        X1 <- predict(gam1, type = "lpmatrix")
        if (!is.null(indexTeq1) && k1.tvc != 0) {
          if (range(X1[, indexTeq1])[1] < 0) 
            stop("Check design matrix for smooth(s) of tvc term(s) in eq. 1.")
        }
        Xd1 <- Xdpred(gam1, data, v1[1])
        gam1$y <- data[, v1[1]]
        st.v1 <- c(gam1$coefficients)
        if (!is.null(indexTeq1)) {
          st.v1[mono.sm.pos1] <- exp(st.v1[mono.sm.pos1])
          while (range(Xd1 %*% st.v1)[1] < 0) st.v1[indexTeq1] <- 0.999 * 
              st.v1[indexTeq1]
          gam1$coefficients <- gam1$coefficients.t <- st.v1
          gam1$coefficients.t[mono.sm.pos1] <- exp(gam1$coefficients.t[mono.sm.pos1])
        }
      }
      gam1$formula <- formula.eq1r
      lsgam1 <- length(gam1$smooth)
      y1 <- y1.test
      if (margins[1] %in% c("LN")) 
        y1 <- log(y1)
      attr(data, "terms") <- NULL
      if (!(surv == TRUE && margins[1] %in% bl)) {
        names(gam1$model)[1] <- as.character(formula.eq1r[2])
        X1 <- predict(gam1, type = "lpmatrix")
        l.sp1 <- length(gam1$sp)
        sp1 <- gam1$sp
      }
      gp1 <- gam1$nsdf
      X1.d2 <- dim(X1)[2]
      form.eq12R <- form.eq12(formula.eq2, data, v2, margins[2], 
                              m1d, m2d)
      formula.eq2 <- form.eq12R$formula.eq1
      formula.eq2r <- form.eq12R$formula.eq1r
      y2 <- form.eq12R$y1
      y2.test <- form.eq12R$y1.test
      y2m <- form.eq12R$y1m
      if (surv == FALSE) 
        gam2 <- eval(substitute(gam(formula.eq2, gamma = infl.fac, 
                                    weights = weights, data = data, knots = knots, 
                                    drop.unused.levels = drop.unused.levels), 
                                list(weights = weights)))
      if (surv == TRUE && !(margins[2] %in% bl)) 
        gam2 <- eval(substitute(gam(formula.eq2, gamma = infl.fac, 
                                    weights = weights * cens2, data = data, knots = knots, 
                                    drop.unused.levels = drop.unused.levels), 
                                list(weights = weights, cens2 = cens2)))
      if (surv == TRUE && margins[2] %in% bl) {
        surv.flex <- TRUE
        f.eq2 <- form.eq12R$f.eq1
        data$urcfcphmwicu <- seq(-10, 10, length.out = dim(data)[1])
        tempb <- eval(substitute(gam(f.eq2, family = cox.ph(), 
                                     data = data, weights = cens2, drop.unused.levels = drop.unused.levels), 
                                 list(cens2 = cens2)))
        data$Sh <- as.vector(mm(predict(tempb, type = "response"), 
                                min.pr = min.pr, max.pr = max.pr))
        cens22 <- ifelse(cens2 == 0, 1e-07, cens2)
        gam2 <- eval(substitute(scam(formula.eq2, gamma = infl.fac, 
                                     weights = weights * cens22, data = data), 
                                list(weights = weights, cens22 = cens22)))
        lsgam2 <- length(gam2$smooth)
        if (lsgam2 == 0) 
          stop("You must use at least a monotonic smooth function of time in the second equation.")
        clsm <- ggr <- NA
        for (i in 1:lsgam2) {
          clsm[i] <- class(gam2$smooth[[i]])[1]
        }
        if (sum(as.numeric(clsm[1] %in% c("mpi.smooth"))) == 
            0) 
          stop("You must have a monotonic smooth of time and it has to be the first to be included.")
        l.sp2 <- length(gam2$sp)
        if (l.sp2 != 0) 
          sp2 <- gam2$sp
        sp2[1] <- 1
        gam.call <- gam2$call
        gam.call$sp <- sp2
        gam2 <- eval(gam.call)
        j <- 1
        for (i in 1:lsgam2) {
          if (max(as.numeric(grepl(v2[1], gam2$smooth[[i]]$term))) != 
              0 && clsm[i] == "mpi.smooth") 
            mono.sm.pos2 <- c(mono.sm.pos2, c(gam2$smooth[[i]]$first.para:gam2$smooth[[i]]$last.para))
        }
        X2 <- predict(gam2, type = "lpmatrix")
        if (!is.null(indexTeq2) && k2.tvc != 0) {
          if (range(X2[, indexTeq2])[1] < 0) 
            stop("Check design matrix for smooth(s) of tvc term(s) in eq. 2.")
        }
        Xd2 <- Xdpred(gam2, data, v2[1])
        gam2$y <- data[, v2[1]]
        st.v2 <- c(gam2$coefficients)
        if (!is.null(indexTeq2)) {
          st.v2[mono.sm.pos2] <- exp(st.v2[mono.sm.pos2])
          while (range(Xd2 %*% st.v2)[1] < 0) st.v2[indexTeq2] <- 0.999 * 
              st.v2[indexTeq2]
          gam2$coefficients <- gam2$coefficients.t <- st.v2
          gam2$coefficients.t[mono.sm.pos2] <- exp(gam2$coefficients.t[mono.sm.pos2])
        }
      }
      gam2$formula <- formula.eq2r
      lsgam2 <- length(gam2$smooth)
      y2 <- y2.test
      if (margins[2] %in% c("LN")) 
        y2 <- log(y2)
      attr(data, "terms") <- NULL
      if (!(surv == TRUE && margins[2] %in% bl)) {
        names(gam2$model)[1] <- as.character(formula.eq2r[2])
        X2 <- predict(gam2, type = "lpmatrix")
        l.sp2 <- length(gam2$sp)
        sp2 <- gam2$sp
      }
      gp2 <- gam2$nsdf
      X2.d2 <- dim(X2)[2]
      res1 <- residuals(gam1)
      res2 <- residuals(gam2)
      ass.s <- cor(res1, res2, method = "kendall")
      ass.s <- sign(ass.s) * ifelse(abs(ass.s) > 0.9, 
                                    0.9, abs(ass.s))
      i.rho <- ass.dp(ass.s, BivD, scc, sccn, nCa)
      dof.st <- log(dof - 2)
      names(dof.st) <- "dof.star"
      if (!(margins[1] %in% c(m1d, bl))) {
        start.snR <- startsn(margins[1], y1)
        log.sig2.1 <- start.snR$log.sig2.1
        names(log.sig2.1) <- "sigma1.star"
        if (margins[1] %in% c(m3)) {
          log.nu.1 <- start.snR$log.nu.1
          names(log.nu.1) <- "nu.1.star"
        }
      }
      if (!(margins[2] %in% c(m1d, bl))) {
        start.snR <- startsn(margins[2], y2)
        log.sig2.2 <- start.snR$log.sig2.1
        names(log.sig2.2) <- "sigma2.star"
        if (margins[2] %in% c(m3)) {
          log.nu.2 <- start.snR$log.nu.1
          names(log.nu.2) <- "nu.2.star"
        }
      }
      vo <- list(gam1 = gam1, gam2 = gam2, i.rho = i.rho, 
                 log.sig2.2 = log.sig2.2, log.nu.2 = log.nu.2, 
                 log.nu.1 = log.nu.1, log.sig2.1 = log.sig2.1, 
                 dof.st = dof.st, n = n, drop.unused.levels = drop.unused.levels)
      start.v <- overall.sv(margins, M, vo)
      if (l.flist > 2) {
        overall.svGR <- overall.svG(formula, data, ngc = 2, 
                                    margins, M, vo, gam1, gam2, knots = knots)
        start.v = overall.svGR$start.v
        X3 = overall.svGR$X3
        X4 = overall.svGR$X4
        X5 = overall.svGR$X5
        X6 = overall.svGR$X6
        X7 = overall.svGR$X7
        X8 = overall.svGR$X8
        X3.d2 = overall.svGR$X3.d2
        X4.d2 = overall.svGR$X4.d2
        X5.d2 = overall.svGR$X5.d2
        X6.d2 = overall.svGR$X6.d2
        X7.d2 = overall.svGR$X7.d2
        X8.d2 = overall.svGR$X8.d2
        gp3 = overall.svGR$gp3
        gp4 = overall.svGR$gp4
        gp5 = overall.svGR$gp5
        gp6 = overall.svGR$gp6
        gp7 = overall.svGR$gp7
        gp8 = overall.svGR$gp8
        gam3 = overall.svGR$gam3
        gam4 = overall.svGR$gam4
        gam5 = overall.svGR$gam5
        gam6 = overall.svGR$gam6
        gam7 = overall.svGR$gam7
        gam8 = overall.svGR$gam8
        l.sp3 = overall.svGR$l.sp3
        l.sp4 = overall.svGR$l.sp4
        l.sp5 = overall.svGR$l.sp5
        l.sp6 = overall.svGR$l.sp6
        l.sp7 = overall.svGR$l.sp7
        l.sp8 = overall.svGR$l.sp8
        sp3 = overall.svGR$sp3
        sp4 = overall.svGR$sp4
        sp5 = overall.svGR$sp5
        sp6 = overall.svGR$sp6
        sp7 = overall.svGR$sp7
        sp8 = overall.svGR$sp8
      }
      GAM <- list(gam1 = gam1, gam2 = gam2, gam3 = gam3, 
                  gam4 = gam4, gam5 = gam5, gam6 = gam6, gam7 = gam7, 
                  gam8 = gam8, gam9 = gam9)
      aregams <- FALSE
      if ((l.sp1 != 0 || l.sp2 != 0 || l.sp3 != 0 || l.sp4 != 
           0 || l.sp5 != 0 || l.sp6 != 0 || l.sp7 != 0 || 
           l.sp8 != 0) && fp == FALSE) {
        L.GAM <- list(l.gam1 = length(gam1$coefficients), 
                      l.gam2 = length(gam2$coefficients), l.gam3 = length(gam3$coefficients), 
                      l.gam4 = length(gam4$coefficients), l.gam5 = length(gam5$coefficients), 
                      l.gam6 = length(gam6$coefficients), l.gam7 = length(gam7$coefficients), 
                      l.gam8 = length(gam8$coefficients), l.gam9 = 0)
        L.SP <- list(l.sp1 = l.sp1, l.sp2 = l.sp2, l.sp3 = l.sp3, 
                     l.sp4 = l.sp4, l.sp5 = l.sp5, l.sp6 = l.sp6, 
                     l.sp7 = l.sp7, l.sp8 = l.sp8, l.sp9 = l.sp9)
        sp <- c(sp1, sp2, sp3, sp4, sp5, sp6, sp7, sp8, 
                sp9)
        qu.mag <- S.m(GAM, L.SP, L.GAM)
        aregams <- TRUE
      }
      if (missing(parscale)) 
        parscale <- 1
      respvec <- respvec2 <- respvec3 <- list(y1 = y1, 
                                              y2 = y2, y1.y2 = NULL, y1.cy2 = NULL, cy1.y2 = NULL, 
                                              cy1.cy2 = NULL, cy1 = NULL, cy = NULL, univ = 0)
      my.env <- new.env()
      my.env$signind <- 1
      lsgam3 <- length(gam3$smooth)
      lsgam4 <- length(gam4$smooth)
      lsgam5 <- length(gam5$smooth)
      lsgam6 <- length(gam6$smooth)
      lsgam7 <- length(gam7$smooth)
      lsgam8 <- length(gam8$smooth)
      lsgam9 <- length(gam9$smooth)
      indUR <- indUL <- indUI <- indUU <- indRR <- indRL <- indRI <- indRU <- indLR <- indLL <- indLI <- indLU <- indIR <- indIL <- indII <- indIU <- rep(0, 
                                                                                                                                                          n)
      if (surv == TRUE && dep.cens == FALSE) {
        if ((surv == TRUE && margins[1] %in% bl && margins[2] %in% 
             bl && !is.factor(cens1) && !is.factor(cens2)) || 
            (surv == TRUE && margins[1] %in% m2 && margins[2] %in% 
             m2)) {
          c11 <- cens1 * cens2
          c10 <- cens1 * (1 - cens2)
          c01 <- (1 - cens1) * cens2
          c00 <- (1 - cens1) * (1 - cens2)
        }
        if (surv == TRUE && margins[1] %in% c(m2, m3) && 
            margins[2] %in% bl) {
          c11 <- cens2
          c10 <- 1 - cens2
          c01 <- NULL
          c00 <- NULL
        }
        if (!is.null(cens1Mix) && !is.null(cens2Mix)) {
          if (surv == TRUE && margins[1] %in% bl && 
              margins[2] %in% bl && is.factor(cens1Mix) && 
              is.factor(cens2Mix)) {
            gamlssfit <- TRUE
            indUR <- as.numeric(cens1Mix == "U") * as.numeric(cens2Mix == 
                                                                "R")
            indUL <- as.numeric(cens1Mix == "U") * as.numeric(cens2Mix == 
                                                                "L")
            indUI <- as.numeric(cens1Mix == "U") * as.numeric(cens2Mix == 
                                                                "I")
            indUU <- as.numeric(cens1Mix == "U") * as.numeric(cens2Mix == 
                                                                "U")
            indRR <- as.numeric(cens1Mix == "R") * as.numeric(cens2Mix == 
                                                                "R")
            indRL <- as.numeric(cens1Mix == "R") * as.numeric(cens2Mix == 
                                                                "L")
            indRI <- as.numeric(cens1Mix == "R") * as.numeric(cens2Mix == 
                                                                "I")
            indRU <- as.numeric(cens1Mix == "R") * as.numeric(cens2Mix == 
                                                                "U")
            indLR <- as.numeric(cens1Mix == "L") * as.numeric(cens2Mix == 
                                                                "R")
            indLL <- as.numeric(cens1Mix == "L") * as.numeric(cens2Mix == 
                                                                "L")
            indLI <- as.numeric(cens1Mix == "L") * as.numeric(cens2Mix == 
                                                                "I")
            indLU <- as.numeric(cens1Mix == "L") * as.numeric(cens2Mix == 
                                                                "U")
            indIR <- as.numeric(cens1Mix == "I") * as.numeric(cens2Mix == 
                                                                "R")
            indIL <- as.numeric(cens1Mix == "I") * as.numeric(cens2Mix == 
                                                                "L")
            indII <- as.numeric(cens1Mix == "I") * as.numeric(cens2Mix == 
                                                                "I")
            indIU <- as.numeric(cens1Mix == "I") * as.numeric(cens2Mix == 
                                                                "U")
          }
        }
      }
      if (surv == TRUE && dep.cens == TRUE) {
        c11 <- NULL
        c10 <- cens1
        c01 <- cens2
        c00 <- cens3
      }
      my.env$k1 <- k1.tvc
      my.env$k2 <- k2.tvc
      VC <- list(lsgam1 = lsgam1, indexTeq1 = indexTeq1, 
                 indexTeq2 = indexTeq2, lsgam2 = lsgam2, Deq1 = Deq1, 
                 pos.pbeq1 = pos.pbeq1, Deq2 = Deq2, pos.pbeq2 = pos.pbeq2, 
                 lsgam3 = lsgam3, robust = FALSE, sp.fixed = NULL, 
                 lsgam4 = lsgam4, Sl.sf = Sl.sf, sp.method = sp.method, 
                 lsgam5 = lsgam5, K1 = NULL, lsgam6 = lsgam6, 
                 lsgam7 = lsgam7, lsgam8 = lsgam8, lsgam9 = lsgam9, 
                 X1 = X1, X2 = X2, X3 = X3, X4 = X4, X5 = X5, 
                 X6 = X6, X7 = X7, X8 = X8, X1.d2 = X1.d2, X2.d2 = X2.d2, 
                 X3.d2 = X3.d2, X4.d2 = X4.d2, X5.d2 = X5.d2, 
                 X6.d2 = X6.d2, X7.d2 = X7.d2, X8.d2 = X8.d2, 
                 gp1 = gp1, gp2 = gp2, gp3 = gp3, gp4 = gp4, 
                 gp5 = gp5, gp6 = gp6, gp7 = gp7, gp8 = gp8, 
                 l.sp1 = l.sp1, l.sp2 = l.sp2, l.sp3 = l.sp3, 
                 l.sp4 = l.sp4, l.sp5 = l.sp5, l.sp6 = l.sp6, 
                 l.sp7 = l.sp7, l.sp8 = l.sp8, l.sp9 = 0, my.env = my.env, 
                 infl.fac = infl.fac, weights = weights, fp = fp, 
                 gamlssfit = gamlssfit, hess = NULL, Model = "CC", 
                 univ.gamls = FALSE, model = model, end = end, 
                 BivD = BivD, nCa = nCa, copula = copula, copula2 = copula2, 
                 nC = nC, gc.l = gc.l, n = n, extra.regI = extra.regI, 
                 parscale = parscale, margins = margins, Cont = "YES", 
                 ccss = "no", m2 = m2, m3 = m3, m1d = m1d, m2d = m2d, 
                 m3d = m3d, bl = bl, triv = FALSE, y1m = y1m, 
                 y2m = y2m, tc = t.c, i.rho = i.rho, dof = dof, 
                 dof.st = dof.st, BivD2 = BivD2, cta = cta, ct = ct, 
                 zerov = -10, c11 = c11, c10 = c10, c01 = c01, 
                 c00 = c00, indUR = indUR, indUL = indUL, indUI = indUI, 
                 indUU = indUU, indRR = indRR, indRL = indRL, 
                 indRI = indRI, indRU = indRU, indLR = indLR, 
                 indLL = indLL, indLI = indLI, indLU = indLU, 
                 indIR = indIR, indIL = indIL, indII = indII, 
                 indIU = indIU, surv = surv, Xd1 = Xd1, Xd2 = Xd2, 
                 mono.sm.pos1 = mono.sm.pos1, mono.sm.pos2 = mono.sm.pos2, 
                 surv.flex = surv.flex, mono.sm.pos = mono.sm.pos, 
                 gp2.inf = NULL, informative = "no", zero.tol = 0.01, 
                 min.dn = min.dn, min.pr = min.pr, max.pr = max.pr)
      if (gc.l == TRUE) 
        gc()
      if (gamlssfit == TRUE) {
        type.cens1 <- type.cens2 <- "R"
        surv1 <- surv2 <- surv
        form.gamlR <- form.gaml(formula, l.flist, M)
        if (surv == TRUE && margins[1] %in% c(m2, m3) && 
            margins[2] %in% bl) 
          surv1 <- FALSE
        if (surv == TRUE && margins[1] %in% bl && margins[2] %in% 
            bl && is.factor(cens1Mix) && is.factor(cens2Mix)) {
          cens1 <- cens1Mix
          cens2 <- cens2Mix
          type.cens1 <- type.cens2 <- "mixed"
          M$type.cens1 = type.cens1
          M$type.cens2 = type.cens2
        }
        gamlss1 <- eval(substitute(gamlss(form.gamlR$formula.gamlss1, 
                                          data = data, weights = weights, subset = subset, 
                                          margin = margins[1], surv = surv1, cens = cens1, 
                                          type.cens = type.cens1, upperB = upperBt1, 
                                          infl.fac = infl.fac, rinit = rinit, rmax = rmax, 
                                          iterlimsp = iterlimsp, tolsp = tolsp, gc.l = gc.l, 
                                          parscale = 1, extra.regI = extra.regI, k.tvc = k1.tvc, 
                                          drop.unused.levels = drop.unused.levels), 
                                   list(weights = weights, cens1 = cens1)))
        gamlss2 <- eval(substitute(gamlss(form.gamlR$formula.gamlss2, 
                                          data = data, weights = weights, subset = subset, 
                                          margin = margins[2], surv = surv2, cens = cens2, 
                                          type.cens = type.cens2, upperB = upperBt2, 
                                          infl.fac = infl.fac, rinit = rinit, rmax = rmax, 
                                          iterlimsp = iterlimsp, tolsp = tolsp, gc.l = gc.l, 
                                          parscale = 1, extra.regI = extra.regI, k.tvc = k2.tvc, 
                                          drop.unused.levels = drop.unused.levels), 
                                   list(weights = weights, cens2 = cens2)))
        SP <- list(sp1 = sp1, sp2 = sp2, sp3 = sp3, 
                   sp4 = sp4, sp5 = sp5, sp6 = sp6, sp7 = sp7, 
                   sp8 = sp8)
        gamls.upsvR <- gamls.upsv(gamlss1, gamlss2, 
                                  margins, M, l.flist, nstv = names(start.v), 
                                  VC, GAM, SP)
        sp <- gamls.upsvR$sp
        start.v <- gamls.upsvR$start.v
        VC$X1 <- gamlss1$VC$X1
        VC$Xd1 <- gamlss1$VC$Xd1
        VC$X1.2 <- gamlss1$VC$X2
        VC$X2 <- gamlss2$VC$X1
        VC$Xd2 <- gamlss2$VC$Xd1
        VC$X2.2 <- gamlss2$VC$X2
        rangeSurv1 <- gamlss1$rangeSurv
        rangeSurv2 <- gamlss2$rangeSurv
      }
      
      
      func.opt <- bdiscrdiscr_ARB
      
      if (aregams) {
        qu.mag$kappad <- TRUE
      }
      

      start.v_use <- c(start.v,
                       logit(0.5*mean(respvec$y1 == 0)+0.01),
                       logit(0.5*mean(respvec$y2 == 0)+0.01))

      SemiParFit <- SemiParBIV.fit(func.opt = func.opt, 
                                   start.v = start.v_use, rinit = rinit, rmax = rmax, 
                                   iterlim = 10000, iterlimsp = iterlimsp, tolsp = tolsp, 
                                   respvec = respvec, VC = VC, sp = sp, qu.mag = qu.mag)
      SemiParFit.p <- copulaReg.fit.post(SemiParFit = SemiParFit, 
                                         VC = VC, GAM)
      y1.m <- y1
      if (margins[1] == "LN") 
        y1.m <- exp(y1)
      y2.m <- y2
      if (margins[2] == "LN") 
        y2.m <- exp(y2)
      SemiParFit <- SemiParFit.p$SemiParFit
      if (gc.l == TRUE) 
        gc()
      

      SemiParFit$fit$e.v <- round(min(eigen(SemiParFit$fit$hessian, symmetric = TRUE, only.values = TRUE)$values), 6)
      SemiParFit$fit$gradi <- round(max(abs(SemiParFit$fit$gradient)), 1)
      
      
      cov.c(SemiParFit)
      
      gam1$call$data <- gam2$call$data <- gam3$call$data <- gam4$call$data <- gam5$call$data <- gam6$call$data <- gam7$call$data <- gam8$call$data <- cl$data
      L <- list(fit = SemiParFit$fit, dataset = NULL, 
                n = n, gamlss1 = gamlss1, gamlss2 = gamlss2, 
                formula = formula, robust = FALSE, edf11 = SemiParFit.p$edf11, 
                surv = surv, gam1 = gam1, gam2 = gam2, gam3 = gam3, 
                gam4 = gam4, gam5 = gam5, gam6 = gam6, gam7 = gam7, 
                gam8 = gam8, coefficients = SemiParFit$fit$argument, 
                coef.t = SemiParFit.p$coef.t, iterlimsp = iterlimsp, 
                weights = weights, cens1 = cens1, cens2 = cens2, 
                cens3 = cens3, sp = SemiParFit.p$sp, iter.sp = SemiParFit$iter.sp, 
                l.sp1 = l.sp1, l.sp2 = l.sp2, l.sp3 = l.sp3, 
                l.sp4 = l.sp4, l.sp5 = l.sp5, l.sp6 = l.sp6, 
                l.sp7 = l.sp7, l.sp8 = l.sp8, bl = bl, l.sp9 = l.sp9, 
                gam9 = gam9, fp = fp, iter.if = SemiParFit$iter.if, 
                iter.inner = SemiParFit$iter.inner, theta = SemiParFit.p$theta, 
                theta.a = SemiParFit.p$theta.a, sigma21 = SemiParFit.p$sigma21, 
                sigma22 = SemiParFit.p$sigma22, sigma21.a = SemiParFit.p$sigma21.a, 
                sigma22.a = SemiParFit.p$sigma22.a, sigma1 = SemiParFit.p$sigma21, 
                sigma2 = SemiParFit.p$sigma22, sigma1.a = SemiParFit.p$sigma21.a, 
                sigma2.a = SemiParFit.p$sigma22.a, nu1 = SemiParFit.p$nu1, 
                nu2 = SemiParFit.p$nu2, nu1.a = SemiParFit.p$nu1.a, 
                nu2.a = SemiParFit.p$nu2.a, dof.a = SemiParFit.p$dof.a, 
                dof = SemiParFit.p$dof, X1 = X1, X2 = X2, X3 = X3, 
                X4 = X4, X5 = X5, X6 = X6, X7 = X7, X8 = X8, 
                X1.d2 = X1.d2, X2.d2 = X2.d2, X3.d2 = X3.d2, 
                X4.d2 = X4.d2, X5.d2 = X5.d2, X6.d2 = X6.d2, 
                X7.d2 = X7.d2, X8.d2 = X8.d2, He = SemiParFit.p$He, 
                HeSh = SemiParFit.p$HeSh, Vb = SemiParFit.p$Vb, 
                Ve = SemiParFit.p$Ve, F = SemiParFit.p$F, F1 = SemiParFit.p$F1, 
                t.edf = SemiParFit.p$t.edf, edf = SemiParFit.p$edf, 
                edf1 = SemiParFit.p$edf1, edf2 = SemiParFit.p$edf2, 
                edf3 = SemiParFit.p$edf3, edf4 = SemiParFit.p$edf4, 
                edf5 = SemiParFit.p$edf5, edf6 = SemiParFit.p$edf6, 
                edf7 = SemiParFit.p$edf7, edf8 = SemiParFit.p$edf8, 
                edf1.1 = SemiParFit.p$edf1.1, edf1.2 = SemiParFit.p$edf1.2, 
                edf1.3 = SemiParFit.p$edf1.3, edf1.4 = SemiParFit.p$edf1.4, 
                edf1.5 = SemiParFit.p$edf1.5, edf1.6 = SemiParFit.p$edf1.6, 
                edf1.7 = SemiParFit.p$edf1.7, edf1.8 = SemiParFit.p$edf1.8, 
                R = SemiParFit.p$R, bs.mgfit = SemiParFit$bs.mgfit, 
                conv.sp = SemiParFit$conv.sp, wor.c = SemiParFit$wor.c, 
                eta1 = SemiParFit$fit$eta1, eta2 = SemiParFit$fit$eta2, 
                etad = SemiParFit$fit$etad, etas1 = SemiParFit$fit$etas1, 
                etas2 = SemiParFit$fit$etas2, y1 = y1.m, y2 = y2.m, 
                BivD = BivD, margins = margins, copula = copula, 
                copula2 = copula2, logLik = SemiParFit.p$logLik, 
                nC = nC, respvec = respvec, hess = TRUE, qu.mag = qu.mag, 
                gp1 = gp1, gp2 = gp2, gp3 = gp3, gp4 = gp4, 
                gp5 = gp5, gp6 = gp6, gp7 = gp7, gp8 = gp8, 
                VC = VC, magpp = SemiParFit$magpp, gamlssfit = gamlssfit, 
                Cont = "YES", tau = SemiParFit.p$tau, tau.a = SemiParFit.p$tau.a, 
                l.flist = l.flist, v1 = v1, v2 = v2, triv = FALSE, 
                univar.gamlss = FALSE, BivD2 = BivD2, call = cl, 
                surv = surv, surv.flex = surv.flex, Vb.t = SemiParFit.p$Vb.t, 
                coef.t = SemiParFit.p$coef.t, Model = "CC", 
                model = model)
      if (BivD %in% BivD2) {
        L$teta1 <- SemiParFit$fit$teta1
        L$teta.ind1 <- SemiParFit$fit$teta.ind1
        L$teta2 <- SemiParFit$fit$teta2
        L$teta.ind2 <- SemiParFit$fit$teta.ind2
        L$Cop1 <- SemiParFit$fit$Cop1
        L$Cop2 <- SemiParFit$fit$Cop2
      }
      class(L) <- c("gjrm", "SemiParBIV")
    }
  }
  L
}







bdiscrdiscr_ARB <- function(params, respvec, VC, ps, AT = FALSE) 
{
  
  # now kappas are at the end
  paramlen <- length(params)
  
  kappa1 <- params[paramlen-1]
  kappa2 <- params[paramlen]
  params_old <- params
  params <- params[-c(paramlen-1, paramlen)]
  
  p1 <- p2 <- pdf1 <- pdf2 <- c.copula.be2 <- c.copula.be1 <- c.copula2.be1be2 <- NA
  eta1 <- VC$X1 %*% params[1:VC$X1.d2]
  eta2 <- VC$X2 %*% params[(VC$X1.d2 + 1):(VC$X1.d2 + VC$X2.d2)]
  etad <- etas1 <- etas2 <- l.ln <- NULL
  if (is.null(VC$X3)) {
    sigma21.st <- etas1 <- params[(VC$X1.d2 + VC$X2.d2 + 
                                     1)]
    sigma22.st <- etas2 <- params[(VC$X1.d2 + VC$X2.d2 + 
                                     2)]
    teta.st <- etad <- params[(VC$X1.d2 + VC$X2.d2 + 3)]
  }
  if (!is.null(VC$X3)) {
    sigma21.st <- etas1 <- VC$X3 %*% params[(VC$X1.d2 + 
                                               VC$X2.d2 + 1):(VC$X1.d2 + VC$X2.d2 + VC$X3.d2)]
    sigma22.st <- etas2 <- VC$X4 %*% params[(VC$X1.d2 + 
                                               VC$X2.d2 + VC$X3.d2 + 1):(VC$X1.d2 + VC$X2.d2 + 
                                                                           VC$X3.d2 + VC$X4.d2)]
    teta.st <- etad <- VC$X5 %*% params[(VC$X1.d2 + VC$X2.d2 + 
                                           VC$X3.d2 + VC$X4.d2 + 1):(VC$X1.d2 + VC$X2.d2 + 
                                                                       VC$X3.d2 + VC$X4.d2 + VC$X5.d2)]
  }
  sstr1 <- esp.tr(sigma21.st, VC$margins[1])
  sstr2 <- esp.tr(sigma22.st, VC$margins[2])
  sigma21.st <- sstr1$vrb.st
  sigma22.st <- sstr2$vrb.st
  sigma21 <- sstr1$vrb
  sigma22 <- sstr2$vrb
  eta1 <- eta.tr(eta1, VC$margins[1])
  eta2 <- eta.tr(eta2, VC$margins[2])
  resT <- teta.tr(VC, teta.st)
  teta.st1 <- teta.st2 <- teta.st <- resT$teta.st
  teta1 <- teta2 <- teta <- resT$teta
  Cop1 <- Cop2 <- VC$BivD
  nC1 <- nC2 <- VC$nC
  teta.ind1 <- as.logical(c(1, 0, round(runif(VC$n - 2))))
  teta.ind2 <- teta.ind1 == FALSE
  if (!(VC$BivD %in% VC$BivD2) && length(teta.st) > 1) {
    teta.st1 <- teta.st[teta.ind1]
    teta.st2 <- teta.st[teta.ind2]
    teta1 <- teta[teta.ind1]
    teta2 <- teta[teta.ind2]
  }
  if (VC$BivD %in% VC$BivD2) {
    if (VC$BivD %in% VC$BivD2[c(1:4, 13:16)]) 
      teta.ind1 <- ifelse(VC$my.env$signind * teta > exp(VC$zerov), 
                          TRUE, FALSE)
    if (VC$BivD %in% VC$BivD2[5:12]) 
      teta.ind1 <- ifelse(VC$my.env$signind * teta > exp(VC$zerov) + 
                            1, TRUE, FALSE)
    teta.ind2 <- teta.ind1 == FALSE
    VC$my.env$signind <- ifelse(teta.ind1 == TRUE, 1, -1)
    teta1 <- teta[teta.ind1]
    teta2 <- -teta[teta.ind2]
    teta.st1 <- teta.st[teta.ind1]
    teta.st2 <- teta.st[teta.ind2]
    if (length(teta) == 1) 
      teta.ind2 <- teta.ind1 <- rep(TRUE, VC$n)
    Cop1Cop2R <- Cop1Cop2(VC$BivD)
    Cop1 <- Cop1Cop2R$Cop1
    Cop2 <- Cop1Cop2R$Cop2
    nC1 <- VC$ct[which(VC$ct[, 1] == Cop1), 2]
    nC2 <- VC$ct[which(VC$ct[, 1] == Cop2), 2]
  }
  dHs1 <- distrHsDiscr(respvec$y1, eta1, sigma21, sigma21.st, 
                       nu = 1, nu.st = 1, margin2 = VC$margins[1], naive = FALSE, 
                       y2m = VC$y1m, min.dn = VC$min.dn, min.pr = VC$min.pr, 
                       max.pr = VC$max.pr)
  dHs2 <- distrHsDiscr(respvec$y2, eta2, sigma22, sigma22.st, 
                       nu = 1, nu.st = 1, margin2 = VC$margins[2], naive = FALSE, 
                       y2m = VC$y2m, min.dn = VC$min.dn, min.pr = VC$min.pr, 
                       max.pr = VC$max.pr)
  pdf1 <- dHs1$pdf2
  pdf2 <- dHs2$pdf2
  p1 <- dHs1$p2
  p2 <- dHs2$p2
  C11 <- C01 <- C10 <- C00 <- NA
  if (length(teta1) != 0) {
    C11[teta.ind1] <- mm(BiCDF(p1[teta.ind1], p2[teta.ind1], 
                               nC1, teta1, VC$dof), min.pr = VC$min.pr, max.pr = VC$max.pr)
    C01[teta.ind1] <- mm(BiCDF(mm(p1[teta.ind1] - pdf1[teta.ind1], 
                                  min.pr = VC$min.pr, max.pr = VC$max.pr), p2[teta.ind1], 
                               nC1, teta1, VC$dof), min.pr = VC$min.pr, max.pr = VC$max.pr)
    C10[teta.ind1] <- mm(BiCDF(p1[teta.ind1], mm(p2[teta.ind1] - 
                                                   pdf2[teta.ind1], min.pr = VC$min.pr, max.pr = VC$max.pr), 
                               nC1, teta1, VC$dof), min.pr = VC$min.pr, max.pr = VC$max.pr)
    C00[teta.ind1] <- mm(BiCDF(mm(p1[teta.ind1] - pdf1[teta.ind1], 
                                  min.pr = VC$min.pr, max.pr = VC$max.pr), mm(p2[teta.ind1] - 
                                                                                pdf2[teta.ind1], min.pr = VC$min.pr, max.pr = VC$max.pr), 
                               nC1, teta1, VC$dof), min.pr = VC$min.pr, max.pr = VC$max.pr)
  }
  if (length(teta2) != 0) {
    C11[teta.ind2] <- mm(BiCDF(p1[teta.ind2], p2[teta.ind2], 
                               nC2, teta2, VC$dof), min.pr = VC$min.pr, max.pr = VC$max.pr)
    C01[teta.ind2] <- mm(BiCDF(mm(p1[teta.ind2] - pdf1[teta.ind2], 
                                  min.pr = VC$min.pr, max.pr = VC$max.pr), p2[teta.ind2], 
                               nC2, teta2, VC$dof), min.pr = VC$min.pr, max.pr = VC$max.pr)
    C10[teta.ind2] <- mm(BiCDF(p1[teta.ind2], mm(p2[teta.ind2] - 
                                                   pdf2[teta.ind2], min.pr = VC$min.pr, max.pr = VC$max.pr), 
                               nC2, teta2, VC$dof), min.pr = VC$min.pr, max.pr = VC$max.pr)
    C00[teta.ind2] <- mm(BiCDF(mm(p1[teta.ind2] - pdf1[teta.ind2], 
                                  min.pr = VC$min.pr, max.pr = VC$max.pr), mm(p2[teta.ind2] - 
                                                                                pdf2[teta.ind2], min.pr = VC$min.pr, max.pr = VC$max.pr), 
                               nC2, teta2, VC$dof), min.pr = VC$min.pr, max.pr = VC$max.pr)
  }
  E <- mm(C11 - C01 - C10 + C00, min.pr = VC$min.pr, max.pr = VC$max.pr)
  l.par <- VC$weights * log(E)
  if (length(teta1) != 0) {
    dHC11F <- copgHs(p1[teta.ind1], p2[teta.ind1], eta1 = NULL, 
                     eta2 = NULL, teta1, teta.st1, Cop1, VC$dof, min.dn = VC$min.dn, 
                     min.pr = VC$min.pr, max.pr = VC$max.pr)
    dHC01F <- copgHs(mm(p1[teta.ind1] - pdf1[teta.ind1], 
                        min.pr = VC$min.pr, max.pr = VC$max.pr), p2[teta.ind1], 
                     eta1 = NULL, eta2 = NULL, teta1, teta.st1, Cop1, 
                     VC$dof, min.dn = VC$min.dn, min.pr = VC$min.pr, 
                     max.pr = VC$max.pr)
    dHC10F <- copgHs(p1[teta.ind1], mm(p2[teta.ind1] - pdf2[teta.ind1], 
                                       min.pr = VC$min.pr, max.pr = VC$max.pr), eta1 = NULL, 
                     eta2 = NULL, teta1, teta.st1, Cop1, VC$dof, min.dn = VC$min.dn, 
                     min.pr = VC$min.pr, max.pr = VC$max.pr)
    dHC00F <- copgHs(mm(p1[teta.ind1] - pdf1[teta.ind1], 
                        min.pr = VC$min.pr, max.pr = VC$max.pr), mm(p2[teta.ind1] - 
                                                                      pdf2[teta.ind1], min.pr = VC$min.pr, max.pr = VC$max.pr), 
                     eta1 = NULL, eta2 = NULL, teta1, teta.st1, Cop1, 
                     VC$dof, min.dn = VC$min.dn, min.pr = VC$min.pr, 
                     max.pr = VC$max.pr)
  }
  if (length(teta2) != 0) {
    dHC11S <- copgHs(p1[teta.ind2], p2[teta.ind2], eta1 = NULL, 
                     eta2 = NULL, teta2, teta.st2, Cop2, VC$dof, min.dn = VC$min.dn, 
                     min.pr = VC$min.pr, max.pr = VC$max.pr)
    dHC01S <- copgHs(mm(p1[teta.ind2] - pdf1[teta.ind2], 
                        min.pr = VC$min.pr, max.pr = VC$max.pr), p2[teta.ind2], 
                     eta1 = NULL, eta2 = NULL, teta2, teta.st2, Cop2, 
                     VC$dof, min.dn = VC$min.dn, min.pr = VC$min.pr, 
                     max.pr = VC$max.pr)
    dHC10S <- copgHs(p1[teta.ind2], mm(p2[teta.ind2] - pdf2[teta.ind2], 
                                       min.pr = VC$min.pr, max.pr = VC$max.pr), eta1 = NULL, 
                     eta2 = NULL, teta2, teta.st2, Cop2, VC$dof, min.dn = VC$min.dn, 
                     min.pr = VC$min.pr, max.pr = VC$max.pr)
    dHC00S <- copgHs(mm(p1[teta.ind2] - pdf1[teta.ind2], 
                        min.pr = VC$min.pr, max.pr = VC$max.pr), mm(p2[teta.ind2] - 
                                                                      pdf2[teta.ind2], min.pr = VC$min.pr, max.pr = VC$max.pr), 
                     eta1 = NULL, eta2 = NULL, teta2, teta.st2, Cop2, 
                     VC$dof, min.dn = VC$min.dn, min.pr = VC$min.pr, 
                     max.pr = VC$max.pr)
  }
  derC11.derp1 <- derC01.derp1 <- derC10.derp1 <- derC00.derp1 <- derC11.derp2 <- derC01.derp2 <- derC10.derp2 <- derC00.derp2 <- derC11.derthet <- derC01.derthet <- derC10.derthet <- derC00.derthet <- derteta.derteta.st <- NA
  if (length(teta1) != 0) {
    derC11.derp1[teta.ind1] <- dHC11F$c.copula.be1
    derC01.derp1[teta.ind1] <- dHC01F$c.copula.be1
    derC10.derp1[teta.ind1] <- dHC10F$c.copula.be1
    derC00.derp1[teta.ind1] <- dHC00F$c.copula.be1
    derC11.derp2[teta.ind1] <- dHC11F$c.copula.be2
    derC01.derp2[teta.ind1] <- dHC01F$c.copula.be2
    derC10.derp2[teta.ind1] <- dHC10F$c.copula.be2
    derC00.derp2[teta.ind1] <- dHC00F$c.copula.be2
    derC11.derthet[teta.ind1] <- dHC11F$c.copula.thet
    derC01.derthet[teta.ind1] <- dHC01F$c.copula.thet
    derC10.derthet[teta.ind1] <- dHC10F$c.copula.thet
    derC00.derthet[teta.ind1] <- dHC00F$c.copula.thet
    derteta.derteta.st[teta.ind1] <- dHC11F$derteta.derteta.st
  }
  if (length(teta2) != 0) {
    derC11.derp1[teta.ind2] <- dHC11S$c.copula.be1
    derC01.derp1[teta.ind2] <- dHC01S$c.copula.be1
    derC10.derp1[teta.ind2] <- dHC10S$c.copula.be1
    derC00.derp1[teta.ind2] <- dHC00S$c.copula.be1
    derC11.derp2[teta.ind2] <- dHC11S$c.copula.be2
    derC01.derp2[teta.ind2] <- dHC01S$c.copula.be2
    derC10.derp2[teta.ind2] <- dHC10S$c.copula.be2
    derC00.derp2[teta.ind2] <- dHC00S$c.copula.be2
    derC11.derthet[teta.ind2] <- dHC11S$c.copula.thet
    derC01.derthet[teta.ind2] <- dHC01S$c.copula.thet
    derC10.derthet[teta.ind2] <- dHC10S$c.copula.thet
    derC00.derthet[teta.ind2] <- dHC00S$c.copula.thet
    derteta.derteta.st[teta.ind2] <- dHC11S$derteta.derteta.st
  }
  derpdf1.dereta1 <- dHs1$derpdf2.dereta2
  derp1.dereta1 <- dHs1$derp2.dereta2
  derp1m1.dereta1 <- derp1.dereta1 - derpdf1.dereta1
  derpdf1.dersigma21.st <- dHs1$derpdf2.dersigma2.st
  derp1.dersigma21.st <- dHs1$derp2.dersigma.st
  derp1m1.dersigma21.st <- derp1.dersigma21.st - derpdf1.dersigma21.st
  derpdf2.dereta2 <- dHs2$derpdf2.dereta2
  derp2.dereta2 <- dHs2$derp2.dereta2
  derp2m1.dereta2 <- derp2.dereta2 - derpdf2.dereta2
  derpdf2.dersigma22.st <- dHs2$derpdf2.dersigma2.st
  derp2.dersigma22.st <- dHs2$derp2.dersigma.st
  derp2m1.dersigma22.st <- derp2.dersigma22.st - derpdf2.dersigma22.st
  fE1 <- (derC11.derp1 - derC10.derp1) * derp1.dereta1 - (derC01.derp1 - 
                                                            derC00.derp1) * derp1m1.dereta1
  fE2 <- (derC11.derp2 - derC01.derp2) * derp2.dereta2 - (derC10.derp2 - 
                                                            derC00.derp2) * derp2m1.dereta2
  fEt <- (derC11.derthet - derC01.derthet - derC10.derthet + 
            derC00.derthet) * derteta.derteta.st
  fE1s <- (derC11.derp1 - derC10.derp1) * derp1.dersigma21.st - 
    (derC01.derp1 - derC00.derp1) * derp1m1.dersigma21.st
  fE2s <- (derC11.derp2 - derC01.derp2) * derp2.dersigma22.st - 
    (derC10.derp2 - derC00.derp2) * derp2m1.dersigma22.st
  dl.dbe1 <- VC$weights * fE1/E
  dl.dbe2 <- VC$weights * fE2/E
  dl.dsigma21.st <- VC$weights * fE1s/E
  dl.dsigma22.st <- VC$weights * fE2s/E
  dl.dteta.st <- VC$weights * fEt/E
  der2C11.derp1p1 <- der2C01.derp1p1 <- der2C10.derp1p1 <- der2C00.derp1p1 <- der2C11.derp2p2 <- der2C01.derp2p2 <- der2C10.derp2p2 <- der2C00.derp2p2 <- der2C11.derp1p2 <- der2C01.derp1p2 <- der2C10.derp1p2 <- der2C00.derp1p2 <- der2C11.derp1t <- der2C01.derp1t <- der2C10.derp1t <- der2C00.derp1t <- der2C11.derp2t <- der2C01.derp2t <- der2C10.derp2t <- der2C00.derp2t <- der2C11.derthet2 <- der2C01.derthet2 <- der2C10.derthet2 <- der2C00.derthet2 <- der2teta.derteta.stteta.st <- NA
  if (length(teta1) != 0) {
    der2C11.derp1p1[teta.ind1] <- dHC11F$c.copula2.be1
    der2C01.derp1p1[teta.ind1] <- dHC01F$c.copula2.be1
    der2C10.derp1p1[teta.ind1] <- dHC10F$c.copula2.be1
    der2C00.derp1p1[teta.ind1] <- dHC00F$c.copula2.be1
    der2C11.derp2p2[teta.ind1] <- dHC11F$c.copula2.be2
    der2C01.derp2p2[teta.ind1] <- dHC01F$c.copula2.be2
    der2C10.derp2p2[teta.ind1] <- dHC10F$c.copula2.be2
    der2C00.derp2p2[teta.ind1] <- dHC00F$c.copula2.be2
    der2C11.derp1p2[teta.ind1] <- dHC11F$c.copula2.be1be2
    der2C01.derp1p2[teta.ind1] <- dHC01F$c.copula2.be1be2
    der2C10.derp1p2[teta.ind1] <- dHC10F$c.copula2.be1be2
    der2C00.derp1p2[teta.ind1] <- dHC00F$c.copula2.be1be2
    der2C11.derp1t[teta.ind1] <- dHC11F$c.copula2.be1t
    der2C01.derp1t[teta.ind1] <- dHC01F$c.copula2.be1t
    der2C10.derp1t[teta.ind1] <- dHC10F$c.copula2.be1t
    der2C00.derp1t[teta.ind1] <- dHC00F$c.copula2.be1t
    der2C11.derp2t[teta.ind1] <- dHC11F$c.copula2.be2t
    der2C01.derp2t[teta.ind1] <- dHC01F$c.copula2.be2t
    der2C10.derp2t[teta.ind1] <- dHC10F$c.copula2.be2t
    der2C00.derp2t[teta.ind1] <- dHC00F$c.copula2.be2t
    der2C11.derthet2[teta.ind1] <- dHC11F$bit1.th2ATE
    der2C01.derthet2[teta.ind1] <- dHC01F$bit1.th2ATE
    der2C10.derthet2[teta.ind1] <- dHC10F$bit1.th2ATE
    der2C00.derthet2[teta.ind1] <- dHC00F$bit1.th2ATE
    der2teta.derteta.stteta.st[teta.ind1] <- dHC11F$der2teta.derteta.stteta.st
  }
  if (length(teta2) != 0) {
    der2C11.derp1p1[teta.ind2] <- dHC11S$c.copula2.be1
    der2C01.derp1p1[teta.ind2] <- dHC01S$c.copula2.be1
    der2C10.derp1p1[teta.ind2] <- dHC10S$c.copula2.be1
    der2C00.derp1p1[teta.ind2] <- dHC00S$c.copula2.be1
    der2C11.derp2p2[teta.ind2] <- dHC11S$c.copula2.be2
    der2C01.derp2p2[teta.ind2] <- dHC01S$c.copula2.be2
    der2C10.derp2p2[teta.ind2] <- dHC10S$c.copula2.be2
    der2C00.derp2p2[teta.ind2] <- dHC00S$c.copula2.be2
    der2C11.derp1p2[teta.ind2] <- dHC11S$c.copula2.be1be2
    der2C01.derp1p2[teta.ind2] <- dHC01S$c.copula2.be1be2
    der2C10.derp1p2[teta.ind2] <- dHC10S$c.copula2.be1be2
    der2C00.derp1p2[teta.ind2] <- dHC00S$c.copula2.be1be2
    der2C11.derp1t[teta.ind2] <- dHC11S$c.copula2.be1t
    der2C01.derp1t[teta.ind2] <- dHC01S$c.copula2.be1t
    der2C10.derp1t[teta.ind2] <- dHC10S$c.copula2.be1t
    der2C00.derp1t[teta.ind2] <- dHC00S$c.copula2.be1t
    der2C11.derp2t[teta.ind2] <- dHC11S$c.copula2.be2t
    der2C01.derp2t[teta.ind2] <- dHC01S$c.copula2.be2t
    der2C10.derp2t[teta.ind2] <- dHC10S$c.copula2.be2t
    der2C00.derp2t[teta.ind2] <- dHC00S$c.copula2.be2t
    der2C11.derthet2[teta.ind2] <- dHC11S$bit1.th2ATE
    der2C01.derthet2[teta.ind2] <- dHC01S$bit1.th2ATE
    der2C10.derthet2[teta.ind2] <- dHC10S$bit1.th2ATE
    der2C00.derthet2[teta.ind2] <- dHC00S$bit1.th2ATE
    der2teta.derteta.stteta.st[teta.ind2] <- dHC11S$der2teta.derteta.stteta.st
  }
  der2pdf1.dereta1 <- dHs1$der2pdf2.dereta2
  der2p1.dereta1eta1 <- dHs1$der2p2.dereta2eta2
  der2p1m1.dereta1eta1 <- der2p1.dereta1eta1 - der2pdf1.dereta1
  der2pdf2.dereta2 <- dHs2$der2pdf2.dereta2
  der2p2.dereta2eta2 <- dHs2$der2p2.dereta2eta2
  der2p2m1.dereta2eta2 <- der2p2.dereta2eta2 - der2pdf2.dereta2
  der2pdf1.dersigma21.st2 <- dHs1$der2pdf2.dersigma2.st2
  der2p1.dersigma21.st2 <- dHs1$der2p2.dersigma2.st2
  der2p1m1.dersigma21.st2 <- der2p1.dersigma21.st2 - der2pdf1.dersigma21.st2
  der2pdf2.dersigma22.st2 <- dHs2$der2pdf2.dersigma2.st2
  der2p2.dersigma22.st2 <- dHs2$der2p2.dersigma2.st2
  der2p2m1.dersigma22.st2 <- der2p2.dersigma22.st2 - der2pdf2.dersigma22.st2
  der2pdf1.dereta1dersigma21.st <- dHs1$der2pdf2.dereta2dersigma2.st
  der2p1.dereta1dersigma21.st <- dHs1$der2p2.dereta2dersigma2.st
  der2p1m1.dereta1dersigma21.st <- der2p1.dereta1dersigma21.st - 
    der2pdf1.dereta1dersigma21.st
  der2pdf2.dereta2dersigma22.st <- dHs2$der2pdf2.dereta2dersigma2.st
  der2p2.dereta2dersigma22.st <- dHs2$der2p2.dereta2dersigma2.st
  der2p2m1.dereta2dersigma22.st <- der2p2.dereta2dersigma22.st - 
    der2pdf2.dereta2dersigma22.st
  d2l.be1.be1 <- -VC$weights * ((E * (der2C11.derp1p1 * derp1.dereta1^2 + 
                                        derC11.derp1 * der2p1.dereta1eta1 - (der2C01.derp1p1 * 
                                                                               derp1m1.dereta1^2 + derC01.derp1 * der2p1m1.dereta1eta1) - 
                                        (der2C10.derp1p1 * derp1.dereta1^2 + derC10.derp1 * 
                                           der2p1.dereta1eta1) + (der2C00.derp1p1 * derp1m1.dereta1^2 + 
                                                                    derC00.derp1 * der2p1m1.dereta1eta1)) - fE1^2)/E^2)
  
  d2l.be2.be2 <- -VC$weights * ((E * (der2C11.derp2p2 * derp2.dereta2^2 + 
                                        derC11.derp2 * der2p2.dereta2eta2 - (der2C01.derp2p2 * 
                                                                               derp2.dereta2^2 + derC01.derp2 * der2p2.dereta2eta2) - 
                                        (der2C10.derp2p2 * derp2m1.dereta2^2 + derC10.derp2 * 
                                           der2p2m1.dereta2eta2) + (der2C00.derp2p2 * derp2m1.dereta2^2 + 
                                                                      derC00.derp2 * der2p2m1.dereta2eta2)) - fE2^2)/E^2)
  d2l.sigma21.sigma21 <- -VC$weights * ((E * (der2C11.derp1p1 * 
                                                derp1.dersigma21.st^2 + derC11.derp1 * der2p1.dersigma21.st2 - 
                                                (der2C01.derp1p1 * derp1m1.dersigma21.st^2 + derC01.derp1 * 
                                                   der2p1m1.dersigma21.st2) - (der2C10.derp1p1 * derp1.dersigma21.st^2 + 
                                                                                 derC10.derp1 * der2p1.dersigma21.st2) + (der2C00.derp1p1 * 
                                                                                                                            derp1m1.dersigma21.st^2 + derC00.derp1 * der2p1m1.dersigma21.st2)) - 
                                           fE1s^2)/E^2)
  d2l.sigma22.sigma22 <- -VC$weights * ((E * (der2C11.derp2p2 * 
                                                derp2.dersigma22.st^2 + derC11.derp2 * der2p2.dersigma22.st2 - 
                                                (der2C01.derp2p2 * derp2.dersigma22.st^2 + derC01.derp2 * 
                                                   der2p2.dersigma22.st2) - (der2C10.derp2p2 * derp2m1.dersigma22.st^2 + 
                                                                               derC10.derp2 * der2p2m1.dersigma22.st2) + (der2C00.derp2p2 * 
                                                                                                                            derp2m1.dersigma22.st^2 + derC00.derp2 * der2p2m1.dersigma22.st2)) - 
                                           fE2s^2)/E^2)
  d2l.rho.rho <- -VC$weights * ((E * (der2C11.derthet2 * derteta.derteta.st^2 + 
                                        derC11.derthet * der2teta.derteta.stteta.st - (der2C01.derthet2 * 
                                                                                         derteta.derteta.st^2 + derC01.derthet * der2teta.derteta.stteta.st) - 
                                        (der2C10.derthet2 * derteta.derteta.st^2 + derC10.derthet * 
                                           der2teta.derteta.stteta.st) + (der2C00.derthet2 * 
                                                                            derteta.derteta.st^2 + derC00.derthet * der2teta.derteta.stteta.st)) - 
                                   fEt^2)/E^2)
  d2l.be1.be2 <- -VC$weights * ((E * (der2C11.derp1p2 * derp1.dereta1 * 
                                        derp2.dereta2 - der2C01.derp1p2 * derp1m1.dereta1 * 
                                        derp2.dereta2 - der2C10.derp1p2 * derp1.dereta1 * derp2m1.dereta2 + 
                                        der2C00.derp1p2 * derp1m1.dereta1 * derp2m1.dereta2) - 
                                   fE1 * fE2)/E^2)
  d2l.be1.sigma22 <- -VC$weights * ((E * (der2C11.derp1p2 * 
                                            derp1.dereta1 * derp2.dersigma22.st - der2C01.derp1p2 * 
                                            derp1m1.dereta1 * derp2.dersigma22.st - der2C10.derp1p2 * 
                                            derp1.dereta1 * derp2m1.dersigma22.st + der2C00.derp1p2 * 
                                            derp1m1.dereta1 * derp2m1.dersigma22.st) - fE1 * fE2s)/E^2)
  d2l.be1.rho <- -VC$weights * ((E * (der2C11.derp1t * derp1.dereta1 * 
                                        derteta.derteta.st - der2C01.derp1t * derp1m1.dereta1 * 
                                        derteta.derteta.st - der2C10.derp1t * derp1.dereta1 * 
                                        derteta.derteta.st + der2C00.derp1t * derp1m1.dereta1 * 
                                        derteta.derteta.st) - fE1 * fEt)/E^2)
  d2l.be1.sigma21 <- -VC$weights * ((E * (der2C11.derp1p1 * 
                                            derp1.dereta1 * derp1.dersigma21.st + derC11.derp1 * 
                                            der2p1.dereta1dersigma21.st - (der2C01.derp1p1 * derp1m1.dereta1 * 
                                                                             derp1m1.dersigma21.st + derC01.derp1 * der2p1m1.dereta1dersigma21.st) - 
                                            (der2C10.derp1p1 * derp1.dereta1 * derp1.dersigma21.st + 
                                               derC10.derp1 * der2p1.dereta1dersigma21.st) + (der2C00.derp1p1 * 
                                                                                                derp1m1.dereta1 * derp1m1.dersigma21.st + derC00.derp1 * 
                                                                                                der2p1m1.dereta1dersigma21.st)) - fE1 * fE1s)/E^2)
  d2l.be2.sigma22 <- -VC$weights * ((E * (der2C11.derp2p2 * 
                                            derp2.dereta2 * derp2.dersigma22.st + derC11.derp2 * 
                                            der2p2.dereta2dersigma22.st - (der2C01.derp2p2 * derp2.dereta2 * 
                                                                             derp2.dersigma22.st + derC01.derp2 * der2p2.dereta2dersigma22.st) - 
                                            (der2C10.derp2p2 * derp2m1.dereta2 * derp2m1.dersigma22.st + 
                                               derC10.derp2 * der2p2m1.dereta2dersigma22.st) + 
                                            (der2C00.derp2p2 * derp2m1.dereta2 * derp2m1.dersigma22.st + 
                                               derC00.derp2 * der2p2m1.dereta2dersigma22.st)) - 
                                       fE2 * fE2s)/E^2)
  d2l.be2.sigma21 <- -VC$weights * ((E * (der2C11.derp1p2 * 
                                            derp1.dersigma21.st * derp2.dereta2 - der2C01.derp1p2 * 
                                            derp1m1.dersigma21.st * derp2.dereta2 - der2C10.derp1p2 * 
                                            derp1.dersigma21.st * derp2m1.dereta2 + der2C00.derp1p2 * 
                                            derp1m1.dersigma21.st * derp2m1.dereta2) - fE1s * fE2)/E^2)
  d2l.be2.rho <- -VC$weights * ((E * (der2C11.derp2t * derp2.dereta2 * 
                                        derteta.derteta.st - der2C01.derp2t * derp2.dereta2 * 
                                        derteta.derteta.st - der2C10.derp2t * derp2m1.dereta2 * 
                                        derteta.derteta.st + der2C00.derp2t * derp2m1.dereta2 * 
                                        derteta.derteta.st) - fE2 * fEt)/E^2)
  d2l.rho.sigma21 <- -VC$weights * ((E * (der2C11.derp1t * 
                                            derp1.dersigma21.st * derteta.derteta.st - der2C01.derp1t * 
                                            derp1m1.dersigma21.st * derteta.derteta.st - der2C10.derp1t * 
                                            derp1.dersigma21.st * derteta.derteta.st + der2C00.derp1t * 
                                            derp1m1.dersigma21.st * derteta.derteta.st) - fE1s * 
                                       fEt)/E^2)
  d2l.rho.sigma22 <- -VC$weights * ((E * (der2C11.derp2t * 
                                            derp2.dersigma22.st * derteta.derteta.st - der2C01.derp2t * 
                                            derp2.dersigma22.st * derteta.derteta.st - der2C10.derp2t * 
                                            derp2m1.dersigma22.st * derteta.derteta.st + der2C00.derp2t * 
                                            derp2m1.dersigma22.st * derteta.derteta.st) - fE2s * 
                                       fEt)/E^2)
  d2l.sigma21.sigma22 <- -VC$weights * ((E * (der2C11.derp1p2 * 
                                                derp1.dersigma21.st * derp2.dersigma22.st - der2C01.derp1p2 * 
                                                derp1m1.dersigma21.st * derp2.dersigma22.st - der2C10.derp1p2 * 
                                                derp1.dersigma21.st * derp2m1.dersigma22.st + der2C00.derp1p2 * 
                                                derp1m1.dersigma21.st * derp2m1.dersigma22.st) - fE1s * 
                                           fE2s)/E^2)
  # if (is.null(VC$X3)) {
  #   G <- -c(colSums(c(dl.dbe1) * VC$X1), colSums(c(dl.dbe2) * 
  #                                                  VC$X2), sum(dl.dsigma21.st), sum(dl.dsigma22.st), 
  #           sum(dl.dteta.st))
  #   be1.be1 <- crossprod(VC$X1 * c(d2l.be1.be1), VC$X1)
  #   be2.be2 <- crossprod(VC$X2 * c(d2l.be2.be2), VC$X2)
  #   be1.be2 <- crossprod(VC$X1 * c(d2l.be1.be2), VC$X2)
  #   be1.rho <- t(t(rowSums(t(VC$X1 * c(d2l.be1.rho)))))
  #   be1.sigma21 <- t(t(rowSums(t(VC$X1 * c(d2l.be1.sigma21)))))
  #   be1.sigma22 <- t(t(rowSums(t(VC$X1 * c(d2l.be1.sigma22)))))
  #   be2.rho <- t(t(rowSums(t(VC$X2 * c(d2l.be2.rho)))))
  #   be2.sigma21 <- t(t(rowSums(t(VC$X2 * c(d2l.be2.sigma21)))))
  #   be2.sigma22 <- t(t(rowSums(t(VC$X2 * c(d2l.be2.sigma22)))))
  #   H <- rbind(cbind(be1.be1, be1.be2, be1.sigma21, be1.sigma22, 
  #                    be1.rho), cbind(t(be1.be2), be2.be2, be2.sigma21, 
  #                                    be2.sigma22, be2.rho), cbind(t(be1.sigma21), t(be2.sigma21), 
  #                                                                 sum(d2l.sigma21.sigma21), sum(d2l.sigma21.sigma22), 
  #                                                                 sum(d2l.rho.sigma21)), cbind(t(be1.sigma22), t(be2.sigma22), 
  #                                                                                              sum(d2l.sigma21.sigma22), sum(d2l.sigma22.sigma22), 
  #                                                                                              sum(d2l.rho.sigma22)), cbind(t(be1.rho), t(be2.rho), 
  #                                                                                                                           sum(d2l.rho.sigma21), sum(d2l.rho.sigma22), sum(d2l.rho.rho)))
  # }
  # if (!is.null(VC$X3)) {
  #   G <- -c(colSums(c(dl.dbe1) * VC$X1), colSums(c(dl.dbe2) * 
  #                                                  VC$X2), colSums(c(dl.dsigma21.st) * VC$X3), colSums(c(dl.dsigma22.st) * 
  #                                                                                                        VC$X4), colSums(c(dl.dteta.st) * VC$X5))
  #   be1.be1 <- crossprod(VC$X1 * c(d2l.be1.be1), VC$X1)
  #   be2.be2 <- crossprod(VC$X2 * c(d2l.be2.be2), VC$X2)
  #   be1.be2 <- crossprod(VC$X1 * c(d2l.be1.be2), VC$X2)
  #   be1.rho <- crossprod(VC$X1 * c(d2l.be1.rho), VC$X5)
  #   be2.rho <- crossprod(VC$X2 * c(d2l.be2.rho), VC$X5)
  #   be1.sigma21 <- crossprod(VC$X1 * c(d2l.be1.sigma21), 
  #                            VC$X3)
  #   be1.sigma22 <- crossprod(VC$X1 * c(d2l.be1.sigma22), 
  #                            VC$X4)
  #   be2.sigma21 <- crossprod(VC$X2 * c(d2l.be2.sigma21), 
  #                            VC$X3)
  #   be2.sigma22 <- crossprod(VC$X2 * c(d2l.be2.sigma22), 
  #                            VC$X4)
  #   sigma21.sigma21 <- crossprod(VC$X3 * c(d2l.sigma21.sigma21), 
  #                                VC$X3)
  #   sigma21.sigma22 <- crossprod(VC$X3 * c(d2l.sigma21.sigma22), 
  #                                VC$X4)
  #   rho.sigma21 <- crossprod(VC$X3 * c(d2l.rho.sigma21), 
  #                            VC$X5)
  #   sigma22.sigma22 <- crossprod(VC$X4 * c(d2l.sigma22.sigma22), 
  #                                VC$X4)
  #   rho.sigma22 <- crossprod(VC$X4 * c(d2l.rho.sigma22), 
  #                            VC$X5)
  #   rho.rho <- crossprod(VC$X5 * c(d2l.rho.rho), VC$X5)
  #   H <- rbind(cbind(be1.be1, be1.be2, be1.sigma21, be1.sigma22, 
  #                    be1.rho), cbind(t(be1.be2), be2.be2, be2.sigma21, 
  #                                    be2.sigma22, be2.rho), cbind(t(be1.sigma21), t(be2.sigma21), 
  #                                                                 sigma21.sigma21, sigma21.sigma22, rho.sigma21), 
  #              cbind(t(be1.sigma22), t(be2.sigma22), t(sigma21.sigma22), 
  #                    sigma22.sigma22, rho.sigma22), cbind(t(be1.rho), 
  #                                                         t(be2.rho), t(rho.sigma21), t(rho.sigma22), 
  #                                                         rho.rho))
  # }
  # res <- -sum(l.par)
  # if (VC$extra.regI == "pC") 
  #   H <- regH(H, type = 1)
  params <- params_old
  
  S.h <- ps$S.h
  if (length(S.h) != 1) {
    S.h1 <- 0.5 * crossprod(params, S.h) %*% params
    S.h2 <- S.h %*% params
  }
  else S.h <- S.h1 <- S.h2 <- 0
  # S.res <- res
  # res <- S.res + S.h1
  # G <- G + S.h2
  # H <- H + S.h
  # if (VC$extra.regI == "sED") 
  #   H <- regH(H, type = 2)
  
  
  ########## Incorporating kappas ##########
  Xd2_len <- VC$X1.d2 + VC$X2.d2 + VC$X3.d2 + VC$X4.d2 + VC$X5.d2
  beta1_idcs <- 1:VC$X1.d2
  beta2_idcs <- (VC$X1.d2+1):(VC$X1.d2 + VC$X2.d2)
  sigma21_idcs <- (VC$X1.d2 + VC$X2.d2 + 1):(VC$X1.d2 + VC$X2.d2 + VC$X3.d2) 
  sigma22_idcs <- (VC$X1.d2 + VC$X2.d2 + VC$X3.d2 + 1):(VC$X1.d2 + VC$X2.d2 + VC$X3.d2 + VC$X4.d2)
  tau_idcs <- (VC$X1.d2 + VC$X2.d2 + VC$X3.d2 + VC$X4.d2 + 1):(Xd2_len)
  
  
  pr1 <- sigmoid(kappa1)
  pr2 <- sigmoid(kappa2)
  
  y1zero <- ifelse(respvec$y1==0, 1, 0)
  y2zero <- ifelse(respvec$y2==0, 1, 0)
  bothzero <- ifelse((respvec$y1==0)&(respvec$y2==0), 1, 0)
  
  
  
  ########## likelihood (length n) vector for copula portion ##########
  dcop <- E  
  
  ##########  likelihood (length n) vector for nb1 marginal portions ##########
  dnb1 <- dHs1$pdf2
  dnb2 <- dHs2$pdf2
  
  ##########  overall likelihood (length n) vector with kappas incorporated ##########
  doverall <- (1-pr1)*(1-pr2)*dcop + (1-pr1)*pr2*dnb1*y2zero + pr1*(1-pr2)*dnb2*y1zero + pr1*pr2*bothzero
  
  
  ########## derivatives of dnb (note it's not derivatives of log likelihood) ##########
  dnb1.dbeta1 <- c(derpdf1.dereta1)*VC$X1
  dnb2.dbeta2 <- c(derpdf2.dereta2)*VC$X2 
  
  dnb1.dsigma21 <- c(derpdf1.dersigma21.st)*VC$X3
  dnb2.dsigma22 <- c(derpdf2.dersigma22.st)*VC$X4
  
  
  dnb1_gradmat <- matrix(0, nrow=nrow(VC$X1), ncol=Xd2_len)
  dnb2_gradmat <- matrix(0, nrow=nrow(VC$X1), ncol=Xd2_len)
  
  dnb1_gradmat[,beta1_idcs] <- dnb1.dbeta1
  dnb1_gradmat[,sigma21_idcs] <- dnb1.dsigma21
  
  dnb2_gradmat[,beta2_idcs] <- dnb2.dbeta2
  dnb2_gradmat[,sigma22_idcs] <- dnb2.dsigma22
  
  
  
  ########## derivatives of dcop (note it's not derivatives of log likelihood) ##########
  
  dcop_gradmat <- matrix(0, nrow=nrow(VC$X1), ncol=Xd2_len)
  
  dcop_gradmat[,beta1_idcs] <- c(fE1)*VC$X1
  dcop_gradmat[,sigma21_idcs] <- c(fE1s)*VC$X3
  dcop_gradmat[,beta2_idcs] <- c(fE2)*VC$X2
  dcop_gradmat[,sigma22_idcs] <- c(fE2s)*VC$X4
  dcop_gradmat[,tau_idcs] <- c(fEt)*VC$X5
  
  
  ########## kappa portion of gradient ##########
  dpr1.dkappa1 <- sigmoid(kappa1)*(1-sigmoid(kappa1))
  
  doverall.dkappa1 <- dpr1.dkappa1*(-(1-pr2)*dcop - pr2*dnb1*y2zero + (1-pr2)*dnb2*y1zero + pr2*bothzero)
  dl.dkappa1 <- c(1/doverall)*doverall.dkappa1
  
  dpr2.dkappa2 <- sigmoid(kappa2)*(1-sigmoid(kappa2))
  
  doverall.dkappa2 <- dpr2.dkappa2*(-(1-pr1)*dcop + (1-pr1)*dnb1*y2zero - pr1*dnb2*y1zero + pr1*bothzero)
  dl.dkappa2 <- c(1/doverall)*doverall.dkappa2
  
  ########## rest of gradient (first grad of lik then of loglik) ##########
  grad_concat <- array(c(dcop_gradmat, c(y2zero)*dnb1_gradmat, c(y1zero)*dnb2_gradmat), dim = c(nrow(dcop_gradmat), ncol(dcop_gradmat), 3))
  
  lik_grad <- apply(X=grad_concat, MARGIN=c(1,2), FUN=function(a, p1, p2) (1-p1)*(1-p2)*a[1] + (1-p1)*p2*a[2] + p1*(1-p2)*a[3], p1=pr1, p2=pr2)
  
  loglik_grad <- c(1/doverall)*lik_grad
  
  
  ########## 2nd derivatives of cop ##########
  twomat_array <- function(mat1, mat2, vec){
    if (nrow(mat1) != nrow(mat2)){
      stop("nrow(mat1) must equal nrow(mat2)")
    }
    if (nrow(mat1) != length(vec)){
      stop("nrow(mat1) must equal length(vec)")
    }
    
    arrdims <- c(ncol(mat1), ncol(mat2), nrow(mat1))
    
    mat1 <- c(vec)*mat1
    
    array(sapply(1:nrow(mat1), function(a) outer(mat1[a,], mat2[a,], "*")), dim=arrdims)
    
  }
  
  
  
  
  mirrormat <- function(m) {
    m[lower.tri(m)] <- t(m)[lower.tri(m)]
    m
  }  
  
  # dcop_hessarr <- array(NA, dim=c(rep(Xd2_len, 2), nrow(VC$X1)))
  # 
  # dcop_hessarr[beta1_idcs, beta1_idcs,] <- twomat_array(mat1=VC$X1, mat2=VC$X1, vec=c(fE1^2/E - d2l.be1.be1*E))
  # dcop_hessarr[beta2_idcs, beta2_idcs,] <- twomat_array(mat1=VC$X2, mat2=VC$X2, vec=c(fE2^2/E - d2l.be2.be2*E))
  # dcop_hessarr[beta1_idcs, beta2_idcs,] <- twomat_array(mat1=VC$X1, mat2=VC$X2, vec=c(fE1*fE2/E - d2l.be1.be2*E))
  # 
  # dcop_hessarr[sigma21_idcs, sigma21_idcs,] <- twomat_array(mat1=VC$X3, mat2=VC$X3, vec=c(fE1s^2/E - d2l.sigma21.sigma21*E))
  # dcop_hessarr[sigma22_idcs, sigma22_idcs,] <- twomat_array(mat1=VC$X4, mat2=VC$X4, vec=c(fE2s^2/E - d2l.sigma22.sigma22*E))
  # dcop_hessarr[sigma21_idcs, sigma22_idcs,] <- twomat_array(mat1=VC$X3, mat2=VC$X4, vec=c(fE1s*fE2s/E - d2l.sigma21.sigma22*E))
  # 
  # dcop_hessarr[tau_idcs, tau_idcs,] <- twomat_array(mat1=VC$X5, mat2=VC$X5, vec=c(fEt^2/E - d2l.rho.rho*E))
  # 
  # dcop_hessarr[beta1_idcs, sigma21_idcs,] <- twomat_array(mat1=VC$X1, mat2=VC$X3, vec=c(fE1*fE1s/E - d2l.be1.sigma21*E))
  # dcop_hessarr[beta1_idcs, sigma22_idcs,] <- twomat_array(mat1=VC$X1, mat2=VC$X4, vec=c(fE1*fE2s/E - d2l.be1.sigma22*E))
  # dcop_hessarr[beta2_idcs, sigma21_idcs,] <- twomat_array(mat1=VC$X2, mat2=VC$X3, vec=c(fE2*fE1s/E - d2l.be2.sigma21*E))
  # dcop_hessarr[beta2_idcs, sigma22_idcs,] <- twomat_array(mat1=VC$X2, mat2=VC$X4, vec=c(fE2*fE2s/E - d2l.be2.sigma22*E))
  # 
  # dcop_hessarr[beta1_idcs, tau_idcs,] <- twomat_array(mat1=VC$X1, mat2=VC$X5, vec=c(fE1*fEt/E - d2l.be1.rho*E))
  # dcop_hessarr[beta2_idcs, tau_idcs,] <- twomat_array(mat1=VC$X2, mat2=VC$X5, vec=c(fE2*fEt/E - d2l.be2.rho*E))
  # 
  # dcop_hessarr[sigma21_idcs, tau_idcs,] <- twomat_array(mat1=VC$X3, mat2=VC$X5, vec=c(fE1s*fEt/E - d2l.rho.sigma21*E))
  # dcop_hessarr[sigma22_idcs, tau_idcs,] <- twomat_array(mat1=VC$X4, mat2=VC$X5, vec=c(fE2s*fEt/E - d2l.rho.sigma22*E))
  # 
  # dcop_hessarr <- array(apply(dcop_hessarr, 3, mirrormat), c(rep(Xd2_len, 2), nrow(VC$X1)))  
  
  dcop_hessum <- matrix(NA, nrow=Xd2_len, ncol=Xd2_len)
  
  
  dcop_hessum[beta1_idcs, beta1_idcs] <- crossprod(VC$X1, VC$X1*c(fE1^2/E - d2l.be1.be1*E)*c(1/doverall))
  dcop_hessum[beta2_idcs, beta2_idcs] <- crossprod(VC$X2, VC$X2*c(fE2^2/E - d2l.be2.be2*E)*c(1/doverall))
  dcop_hessum[beta1_idcs, beta2_idcs] <- crossprod(VC$X1, VC$X2*c(fE1*fE2/E - d2l.be1.be2*E)*c(1/doverall))
  
  dcop_hessum[sigma21_idcs, sigma21_idcs] <- crossprod(VC$X3, VC$X3*c(fE1s^2/E - d2l.sigma21.sigma21*E)*c(1/doverall))
  dcop_hessum[sigma22_idcs, sigma22_idcs] <- crossprod(VC$X4, VC$X4*c(fE2s^2/E - d2l.sigma22.sigma22*E)*c(1/doverall))
  dcop_hessum[sigma21_idcs, sigma22_idcs] <- crossprod(VC$X3, VC$X4*c(fE1s*fE2s/E - d2l.sigma21.sigma22*E)*c(1/doverall))
  
  dcop_hessum[tau_idcs, tau_idcs] <- crossprod(VC$X5, VC$X5*c(fEt^2/E - d2l.rho.rho*E)*c(1/doverall))
  
  dcop_hessum[beta1_idcs, sigma21_idcs] <- crossprod(VC$X1, VC$X3*c(fE1*fE1s/E - d2l.be1.sigma21*E)*c(1/doverall))
  dcop_hessum[beta1_idcs, sigma22_idcs] <- crossprod(VC$X1, VC$X4*c(fE1*fE2s/E - d2l.be1.sigma22*E)*c(1/doverall))
  dcop_hessum[beta2_idcs, sigma21_idcs] <- crossprod(VC$X2, VC$X3*c(fE2*fE1s/E - d2l.be2.sigma21*E)*c(1/doverall))
  dcop_hessum[beta2_idcs, sigma22_idcs] <- crossprod(VC$X2, VC$X4*c(fE2*fE2s/E - d2l.be2.sigma22*E)*c(1/doverall))
  
  dcop_hessum[beta1_idcs, tau_idcs] <- crossprod(VC$X1, VC$X5*c(fE1*fEt/E - d2l.be1.rho*E)*c(1/doverall))
  dcop_hessum[beta2_idcs, tau_idcs] <- crossprod(VC$X2, VC$X5*c(fE2*fEt/E - d2l.be2.rho*E)*c(1/doverall))
  
  dcop_hessum[sigma21_idcs, tau_idcs] <- crossprod(VC$X3, VC$X5*c(fE1s*fEt/E - d2l.rho.sigma21*E)*c(1/doverall))
  dcop_hessum[sigma22_idcs, tau_idcs] <- crossprod(VC$X4, VC$X5*c(fE2s*fEt/E - d2l.rho.sigma22*E)*c(1/doverall))
  
  
  
  dcop_hessum <- mirrormat(dcop_hessum)
  
  
  
  ########## 2nd derivatives of nb ##########
  
  # dnb1_hessarr[beta1_idcs, beta1_idcs,] <- twomat_array(mat1=VC$X1, mat2=VC$X1, vec=c(der2pdf1.dereta1))
  # dnb1_hessarr[sigma21_idcs, sigma21_idcs,] <- twomat_array(mat1=VC$X3, mat2=VC$X3, vec=c(der2pdf1.dersigma21.st2))
  # dnb1_hessarr[beta1_idcs, sigma21_idcs,] <- twomat_array(mat1=VC$X1, mat2=VC$X3, vec=c(der2pdf1.dereta1dersigma21.st))
  # 
  # dnb2_hessarr[beta2_idcs, beta2_idcs,] <- twomat_array(mat1=VC$X2, mat2=VC$X2, vec=c(der2pdf2.dereta2))
  # dnb2_hessarr[sigma22_idcs, sigma22_idcs,] <- twomat_array(mat1=VC$X4, mat2=VC$X4, vec=c(der2pdf2.dersigma22.st2))
  # dnb2_hessarr[beta2_idcs, sigma22_idcs,] <- twomat_array(mat1=VC$X2, mat2=VC$X4, vec=c(der2pdf2.dereta2dersigma22.st))
  #   
  # dnb1_hessarr <- array(apply(dnb1_hessarr, 3, mirrormat), dim(dcop_hessarr))
  # dnb2_hessarr <- array(apply(dnb2_hessarr, 3, mirrormat), dim(dcop_hessarr))  
  # dnb1_hessarr <- array(0, dim=dim(dcop_hessarr))
  # dnb2_hessarr <- array(0, dim=dim(dcop_hessarr))
  
  
  dnb1_hessum <- matrix(0, nrow=nrow(dcop_hessum), ncol=ncol(dcop_hessum))
  dnb2_hessum <- matrix(0, nrow=nrow(dcop_hessum), ncol=ncol(dcop_hessum))
  
  dnb1_hessum[beta1_idcs, beta1_idcs] <- crossprod(VC$X1, VC$X1*c(der2pdf1.dereta1)*c(1/doverall)*y2zero)
  dnb1_hessum[sigma21_idcs, sigma21_idcs] <- crossprod(VC$X3, VC$X3*c(der2pdf1.dersigma21.st2)*c(1/doverall)*y2zero)
  dnb1_hessum[beta1_idcs, sigma21_idcs] <- crossprod(VC$X1, VC$X3*c(der2pdf1.dereta1dersigma21.st)*c(1/doverall)*y2zero)
  
  dnb2_hessum[beta2_idcs, beta2_idcs] <- crossprod(VC$X2, VC$X2*c(der2pdf2.dereta2)*c(1/doverall)*y1zero)
  dnb2_hessum[sigma22_idcs, sigma22_idcs] <- crossprod(VC$X4, VC$X4*c(der2pdf2.dersigma22.st2)*c(1/doverall)*y1zero)
  dnb2_hessum[beta2_idcs, sigma22_idcs] <- crossprod(VC$X2, VC$X4*c(der2pdf2.dereta2dersigma22.st)*c(1/doverall)*y1zero)
  
  dnb1_hessum <- mirrormat(dnb1_hessum)
  dnb2_hessum <- mirrormat(dnb2_hessum)
  
  
  ########## combining the nb1, nb2, and cop hessians to give 2nd der of lik (not doing kappas yet) ##########
  # zero out dnb11_hessarr and dnb2_hessarr according to zinf
  
  # dnb1_hessarr[,,which(y2zero == 0)] <- 0
  # dnb2_hessarr[,,which(y1zero == 0)] <- 0
  # 
  # 
  # loglik_hess <- array(NA, dim(dnb1_hessarr))
  # 
  # for (i in 1:Xd2_len){
  #   for (j in i:Xd2_len){
  #     loglik_hess[i,j,] <- c(-1/doverall^2)*lik_grad[,i]*lik_grad[,j] + c(1/doverall)*((1-pr1)*(1-pr2)*dcop_hessarr[i,j,] + (1-pr1)*pr2*dnb1_hessarr[i,j,] + pr1*(1-pr2)*dnb2_hessarr[i,j,])
  #     loglik_hess[j,i,] <- loglik_hess[i,j,]
  #   }
  # }
  
  # hessum <- apply(loglik_hess, c(1,2), sum)
  
  piece1 <- crossprod(lik_grad, lik_grad*c(-1/doverall^2))
  
  piece2 <- (1-pr1)*(1-pr2)*dcop_hessum + (1-pr1)*pr2*dnb1_hessum + pr1*(1-pr2)*dnb2_hessum
  
  
  
  hessum <- piece1 + piece2
  
  
  ########## 2nd derivatives pertaining to kappa ##########
  ########## these are derivatives of logliks ##########
  dpr1.dkappa1kappa1 <- -exp(kappa1)*(exp(kappa1)-1)/(exp(kappa1)+1)^3
  
  dl.dkappa1kappa1 <- (-1/doverall^2)*doverall.dkappa1^2 + c(1/doverall)*dpr1.dkappa1kappa1*doverall.dkappa1/dpr1.dkappa1
  
  dpr2.dkappa2kappa2 <- -exp(kappa2)*(exp(kappa2)-1)/(exp(kappa2)+1)^3
  
  dl.dkappa2kappa2 <- (-1/doverall^2)*doverall.dkappa2^2 + c(1/doverall)*dpr2.dkappa2kappa2*doverall.dkappa2/dpr2.dkappa2
  
  dl.dkappa1kappa2 <- (-1/doverall^2)*doverall.dkappa1*doverall.dkappa2 + c(1/doverall)*(dcop - dnb1*y2zero - dnb2*y1zero + bothzero)*dpr1.dkappa1*dpr2.dkappa2
  
  
  dl.dkappa1else <- c(-doverall.dkappa1/doverall^2)*lik_grad + 
    c(1/doverall)*apply(X=grad_concat, MARGIN=c(1,2), FUN=function(a, p1, p2) -(1-p2)*a[1] - p2*a[2] + (1-p2)*a[3], p1=pr1, p2=pr2)*dpr1.dkappa1
  
  dl.dkappa2else <- c(-doverall.dkappa2/doverall^2)*lik_grad + 
    c(1/doverall)*apply(X=grad_concat, MARGIN=c(1,2), FUN=function(a, p1, p2) -(1-p1)*a[1] + (1-p1)*a[2] - p1*a[3], p1=pr1, p2=pr2)*dpr2.dkappa2
  
  kappa_2nders <- cbind(colSums(cbind(dl.dkappa1else, dl.dkappa1kappa1, dl.dkappa1kappa2)), colSums(cbind(dl.dkappa2else, dl.dkappa1kappa2, dl.dkappa2kappa2)))
  
  
  ########## form the value, gradient and hessian ##########
  G_temp <- -colSums(cbind(loglik_grad, dl.dkappa1, dl.dkappa2))
  H_temp <- matrix(NA, nrow=Xd2_len+2, ncol=Xd2_len+2)
  H_temp[c(Xd2_len+1, Xd2_len+2),] <- t(kappa_2nders)
  H_temp[,c(Xd2_len+1, Xd2_len+2)] <- kappa_2nders
  H_temp[1:Xd2_len,1:Xd2_len] <- hessum
  H_temp <- -H_temp
  
  l.par_temp <- c(log(doverall))
  res_temp <- -sum(l.par_temp)
  S.res_temp <- res_temp
  
  res_temp <- S.res_temp + S.h1
  if (length(S.h2)>1){
    G_temp <- G_temp + S.h2    
  }
  if (length(S.h)>1){
    H_temp <- H_temp + S.h
  }
  
  print(res_temp)
  
  list(value = res_temp, gradient = G_temp, hessian = H_temp, S.h = S.h, 
       S.h1 = S.h1, S.h2 = S.h2, l = S.res_temp, l.ln = l.ln, l.par = l.par_temp, 
       ps = ps, eta1 = eta1, eta2 = eta2, etad = etad, etas1 = etas1, 
       etas2 = etas2, BivD = VC$BivD)
}








pen <- function (qu.mag, sp, VC, univ, l.splist) 
{
  ma1 <- matrix(0, VC$gp1, VC$gp1)
  if (l.splist$l.sp1 == 0) 
    EQ1P <- adiag(ma1)
  if (l.splist$l.sp1 != 0) {
    ind <- 1:l.splist$l.sp1
    offtemp <- as.numeric(as.factor(qu.mag$off[ind]))
    S1 <- mapply("*", qu.mag$Ss[ind], sp[ind], SIMPLIFY = FALSE)
    if (length(unique(offtemp)) != length(offtemp)) 
      S1 <- SS(offtemp, S1)
    S1 <- do.call(adiag, lapply(S1, unlist))
    EQ1P <- adiag(ma1, S1)
  }
  if (is.null(VC$gp2.inf)) 
    ma2 <- matrix(0, VC$gp2, VC$gp2)
  else ma2 <- matrix(0, VC$gp2.inf, VC$gp2.inf)
  if (l.splist$l.sp2 == 0) 
    EQ2P <- adiag(ma2)
  if (l.splist$l.sp2 != 0) {
    ind <- (l.splist$l.sp1 + 1):(l.splist$l.sp1 + l.splist$l.sp2)
    offtemp <- as.numeric(as.factor(qu.mag$off[ind]))
    S2 <- mapply("*", qu.mag$Ss[ind], sp[ind], SIMPLIFY = FALSE)
    if (length(unique(offtemp)) != length(offtemp)) 
      S2 <- SS(offtemp, S2)
    S2 <- do.call(adiag, lapply(S2, unlist))
    EQ2P <- adiag(ma2, S2)
  }
  if (!is.null(VC$gp3)) {
    EQ4P <- EQ5P <- EQ6P <- EQ7P <- EQ8P <- EQ9P <- NULL
    ma3 <- matrix(0, VC$gp3, VC$gp3)
    if (l.splist$l.sp3 == 0) 
      EQ3P <- adiag(ma3)
    if (l.splist$l.sp3 != 0) {
      ind <- (l.splist$l.sp1 + l.splist$l.sp2 + 1):(l.splist$l.sp1 + 
                                                      l.splist$l.sp2 + l.splist$l.sp3)
      offtemp <- as.numeric(as.factor(qu.mag$off[ind]))
      S3 <- mapply("*", qu.mag$Ss[ind], sp[ind], SIMPLIFY = FALSE)
      if (length(unique(offtemp)) != length(offtemp)) 
        S3 <- SS(offtemp, S3)
      S3 <- do.call(adiag, lapply(S3, unlist))
      EQ3P <- adiag(ma3, S3)
    }
    if (!is.null(VC$gp4)) {
      ma4 <- matrix(0, VC$gp4, VC$gp4)
      if (l.splist$l.sp4 == 0) 
        EQ4P <- adiag(ma4)
      if (l.splist$l.sp4 != 0) {
        ind <- (l.splist$l.sp1 + l.splist$l.sp2 + l.splist$l.sp3 + 
                  1):(l.splist$l.sp1 + l.splist$l.sp2 + l.splist$l.sp3 + 
                        l.splist$l.sp4)
        offtemp <- as.numeric(as.factor(qu.mag$off[ind]))
        S4 <- mapply("*", qu.mag$Ss[ind], sp[ind], SIMPLIFY = FALSE)
        if (length(unique(offtemp)) != length(offtemp)) 
          S4 <- SS(offtemp, S4)
        S4 <- do.call(adiag, lapply(S4, unlist))
        EQ4P <- adiag(ma4, S4)
      }
    }
    if (!is.null(VC$gp5)) {
      ma5 <- matrix(0, VC$gp5, VC$gp5)
      if (l.splist$l.sp5 == 0) 
        EQ5P <- adiag(ma5)
      if (l.splist$l.sp5 != 0) {
        ind <- (l.splist$l.sp1 + l.splist$l.sp2 + l.splist$l.sp3 + 
                  l.splist$l.sp4 + 1):(l.splist$l.sp1 + l.splist$l.sp2 + 
                                         l.splist$l.sp3 + l.splist$l.sp4 + l.splist$l.sp5)
        offtemp <- as.numeric(as.factor(qu.mag$off[ind]))
        S5 <- mapply("*", qu.mag$Ss[ind], sp[ind], SIMPLIFY = FALSE)
        if (length(unique(offtemp)) != length(offtemp)) 
          S5 <- SS(offtemp, S5)
        S5 <- do.call(adiag, lapply(S5, unlist))
        EQ5P <- adiag(ma5, S5)
      }
    }
    if (!is.null(VC$gp6)) {
      ma6 <- matrix(0, VC$gp6, VC$gp6)
      if (l.splist$l.sp6 == 0) 
        EQ6P <- adiag(ma6)
      if (l.splist$l.sp6 != 0) {
        ind <- (l.splist$l.sp1 + l.splist$l.sp2 + l.splist$l.sp3 + 
                  l.splist$l.sp4 + l.splist$l.sp5 + 1):(l.splist$l.sp1 + 
                                                          l.splist$l.sp2 + l.splist$l.sp3 + l.splist$l.sp4 + 
                                                          l.splist$l.sp5 + l.splist$l.sp6)
        offtemp <- as.numeric(as.factor(qu.mag$off[ind]))
        S6 <- mapply("*", qu.mag$Ss[ind], sp[ind], SIMPLIFY = FALSE)
        if (length(unique(offtemp)) != length(offtemp)) 
          S6 <- SS(offtemp, S6)
        S6 <- do.call(adiag, lapply(S6, unlist))
        EQ6P <- adiag(ma6, S6)
      }
    }
    if (!is.null(VC$gp7)) {
      ma7 <- matrix(0, VC$gp7, VC$gp7)
      if (l.splist$l.sp7 == 0) 
        EQ7P <- adiag(ma7)
      if (l.splist$l.sp7 != 0) {
        ind <- (l.splist$l.sp1 + l.splist$l.sp2 + l.splist$l.sp3 + 
                  l.splist$l.sp4 + l.splist$l.sp5 + l.splist$l.sp6 + 
                  1):(l.splist$l.sp1 + l.splist$l.sp2 + l.splist$l.sp3 + 
                        l.splist$l.sp4 + l.splist$l.sp5 + l.splist$l.sp6 + 
                        l.splist$l.sp7)
        offtemp <- as.numeric(as.factor(qu.mag$off[ind]))
        S7 <- mapply("*", qu.mag$Ss[ind], sp[ind], SIMPLIFY = FALSE)
        if (length(unique(offtemp)) != length(offtemp)) 
          S7 <- SS(offtemp, S7)
        S7 <- do.call(adiag, lapply(S7, unlist))
        EQ7P <- adiag(ma7, S7)
      }
    }
    if (!is.null(VC$gp8)) {
      ma8 <- matrix(0, VC$gp8, VC$gp8)
      if (l.splist$l.sp8 == 0) 
        EQ8P <- adiag(ma8)
      if (l.splist$l.sp8 != 0) {
        ind <- (l.splist$l.sp1 + l.splist$l.sp2 + l.splist$l.sp3 + 
                  l.splist$l.sp4 + l.splist$l.sp5 + l.splist$l.sp6 + 
                  l.splist$l.sp7 + 1):(l.splist$l.sp1 + l.splist$l.sp2 + 
                                         l.splist$l.sp3 + l.splist$l.sp4 + l.splist$l.sp5 + 
                                         l.splist$l.sp6 + l.splist$l.sp7 + l.splist$l.sp8)
        offtemp <- as.numeric(as.factor(qu.mag$off[ind]))
        S8 <- mapply("*", qu.mag$Ss[ind], sp[ind], SIMPLIFY = FALSE)
        if (length(unique(offtemp)) != length(offtemp)) 
          S8 <- SS(offtemp, S8)
        S8 <- do.call(adiag, lapply(S8, unlist))
        EQ8P <- adiag(ma8, S8)
      }
    }
    if (!is.null(VC$gp9)) {
      ma9 <- matrix(0, VC$gp9, VC$gp9)
      if (l.splist$l.sp9 == 0) 
        EQ9P <- adiag(ma9)
      if (l.splist$l.sp9 != 0) {
        ind <- (l.splist$l.sp1 + l.splist$l.sp2 + l.splist$l.sp3 + 
                  l.splist$l.sp4 + l.splist$l.sp5 + l.splist$l.sp6 + 
                  l.splist$l.sp7 + l.splist$l.sp8 + 1):(l.splist$l.sp1 + 
                                                          l.splist$l.sp2 + l.splist$l.sp3 + l.splist$l.sp4 + 
                                                          l.splist$l.sp5 + l.splist$l.sp6 + l.splist$l.sp7 + 
                                                          l.splist$l.sp8 + l.splist$l.sp9)
        offtemp <- as.numeric(as.factor(qu.mag$off[ind]))
        S9 <- mapply("*", qu.mag$Ss[ind], sp[ind], SIMPLIFY = FALSE)
        if (length(unique(offtemp)) != length(offtemp)) 
          S9 <- SS(offtemp, S9)
        S9 <- do.call(adiag, lapply(S9, unlist))
        EQ9P <- adiag(ma9, S9)
      }
    }
  }
  else {
    if (VC$Model != "ROY") {
      if (VC$univ.gamls == FALSE) {
        if (VC$margins[1] %in% c(VC$m2, VC$m3) && VC$margins[2] %in% 
            c(VC$m2, VC$m3) && VC$BivD == "T") {
          if (VC$margins[1] %in% VC$m2 && VC$margins[2] %in% 
              VC$m2) {
            EQ3P <- 0
            EQ4P <- 0
            EQ5P <- 0
            EQ6P <- 0
            EQ7P <- EQ8P <- NULL
          }
          if (VC$margins[1] %in% VC$m3 && VC$margins[2] %in% 
              VC$m3) {
            EQ3P <- 0
            EQ4P <- 0
            EQ5P <- 0
            EQ6P <- 0
            EQ7P <- EQ8P <- 0
          }
          if (VC$margins[1] %in% VC$m2 && VC$margins[2] %in% 
              VC$m3) {
            EQ3P <- 0
            EQ4P <- 0
            EQ5P <- 0
            EQ6P <- 0
            EQ7P <- 0
            EQ8P <- NULL
          }
          if (VC$margins[1] %in% VC$m3 && VC$margins[2] %in% 
              VC$m2) {
            EQ3P <- 0
            EQ4P <- 0
            EQ5P <- 0
            EQ6P <- 0
            EQ7P <- 0
            EQ8P <- NULL
          }
        }
        else {
          if (VC$margins[1] %in% c(VC$bl) && VC$Model != 
              "BPO0") {
            EQ3P <- 0
            EQ4P <- NULL
            EQ5P <- NULL
            EQ6P <- EQ7P <- EQ8P <- NULL
          }
          if (VC$margins[1] %in% c(VC$bl, VC$m1d) && 
              VC$margins[2] %in% c(VC$bl, VC$m1d)) {
            EQ3P <- 0
            EQ4P <- NULL
            EQ5P <- NULL
            EQ6P <- EQ7P <- EQ8P <- NULL
          }
          if (VC$margins[1] %in% c(VC$bl, VC$m1d) && 
              VC$margins[2] %in% VC$m2d) {
            EQ3P <- 0
            EQ4P <- 0
            EQ5P <- NULL
            EQ6P <- EQ7P <- EQ8P <- NULL
          }
          if (VC$margins[1] %in% c(VC$bl, VC$m1d) && 
              VC$margins[2] %in% VC$m2) {
            EQ3P <- 0
            EQ4P <- 0
            EQ5P <- NULL
            EQ6P <- EQ7P <- EQ8P <- NULL
          }
          if (VC$margins[1] %in% c(VC$bl, VC$m1d) && 
              VC$margins[2] %in% VC$m3) {
            EQ3P <- 0
            EQ4P <- 0
            EQ5P <- 0
            EQ6P <- EQ7P <- EQ8P <- NULL
          }
          if (VC$margins[1] %in% c(VC$m2, VC$m2d) && 
              VC$margins[2] %in% c(VC$m2, VC$m2d)) {
            EQ3P <- 0
            EQ4P <- 0
            EQ5P <- 0
            EQ6P <- NULL
            EQ7P <- EQ8P <- NULL
          }
          if (VC$margins[1] %in% c(VC$m2, VC$m2d) && 
              VC$margins[2] %in% c(VC$bl)) {
            EQ3P <- 0
            EQ4P <- 0
            EQ5P <- NULL
            EQ6P <- NULL
            EQ7P <- EQ8P <- NULL
          }
          if (VC$margins[1] %in% c(VC$m3) && VC$margins[2] %in% 
              c(VC$bl)) {
            EQ3P <- 0
            EQ4P <- 0
            EQ5P <- 0
            EQ6P <- NULL
            EQ7P <- EQ8P <- NULL
          }
          if (VC$margins[1] %in% VC$m3 && VC$margins[2] %in% 
              VC$m3) {
            EQ3P <- 0
            EQ4P <- 0
            EQ5P <- 0
            EQ6P <- 0
            EQ7P <- 0
            EQ8P <- NULL
          }
          if (VC$margins[1] %in% c(VC$m2, VC$m2d) && 
              VC$margins[2] %in% VC$m3) {
            EQ3P <- 0
            EQ4P <- 0
            EQ5P <- 0
            EQ6P <- 0
            EQ7P <- EQ8P <- NULL
          }
          if (VC$margins[1] %in% VC$m3 && VC$margins[2] %in% 
              c(VC$m2, VC$m2d)) {
            EQ3P <- 0
            EQ4P <- 0
            EQ5P <- 0
            EQ6P <- 0
            EQ7P <- EQ8P <- NULL
          }
          if (VC$Model == "B" && !is.null(VC$theta.fx)) {
            EQ3P <- EQ4P <- EQ5P <- EQ6P <- EQ7P <- EQ8P <- NULL
          }
          if (VC$Model == "BPO0") {
            EQ3P <- EQ4P <- EQ5P <- EQ6P <- EQ7P <- EQ8P <- NULL
          }
        }
      }
    }
  }
  if (VC$triv == FALSE && VC$Model != "ROY") {
    if (univ == 0) 
      S.h <- adiag(EQ1P, EQ2P, EQ3P, EQ4P, EQ5P, EQ6P, 
                   EQ7P, EQ8P)
    if (univ == 2) {
      if (VC$margins[1] %in% c(VC$m1d, VC$bl)) 
        S.h <- adiag(EQ1P)
      if (VC$margins[1] %in% c(VC$bl) && !is.null(VC$gp2.inf)) 
        S.h <- adiag(EQ1P, EQ2P)
      if (VC$margins[1] %in% c(VC$m2, VC$m2d)) 
        S.h <- adiag(EQ1P, EQ2P)
      if (VC$margins[1] %in% VC$m3) 
        S.h <- adiag(EQ1P, EQ2P, EQ3P)
    }
  }
  if (VC$triv == TRUE) {
    S.h <- adiag(EQ1P, EQ2P, EQ3P)
    if (VC$penCor %in% c("unpen") && VC$l.flist == 3) 
      S.h <- adiag(S.h, matrix(0, 3, 3))
    if (VC$penCor %in% c("unpen") && VC$l.flist == 6) 
      S.h <- adiag(EQ1P, EQ2P, EQ3P, EQ4P, EQ5P, EQ6P)
    if (VC$penCor %in% c("ridge")) {
      A <- diag(c(1, 1, 1))
      if (VC$l.sp1 == 0 && VC$l.sp2 == 0 && VC$l.sp3 == 
          0) 
        qu.mag$Ss[[1]] <- A
      if (VC$l.sp1 != 0 || VC$l.sp2 != 0 || VC$l.sp3 != 
          0) 
        qu.mag$Ss[[length(qu.mag$Ss) + 1]] <- A
      S.h <- adiag(S.h, sp[length(sp)] * qu.mag$Ss[[length(qu.mag$Ss)]])
    }
  }
  if (VC$Model == "ROY") {
    if (VC$l.flist == 3) {
      if (VC$margins[2] %in% c(VC$bl) && VC$margins[3] %in% 
          c(VC$bl)) 
        EQ4P <- EQ5P <- 0
      if (VC$margins[2] %in% c(VC$m1d) && VC$margins[3] %in% 
          c(VC$m1d)) 
        EQ4P <- EQ5P <- 0
      if (VC$margins[2] %in% c(VC$m2d) && VC$margins[3] %in% 
          c(VC$m1d)) 
        EQ4P <- EQ5P <- EQ6P <- 0
      if (VC$margins[2] %in% c(VC$m1d) && VC$margins[3] %in% 
          c(VC$m2d)) 
        EQ4P <- EQ5P <- EQ6P <- 0
      if (VC$margins[2] %in% c(VC$m3) && VC$margins[3] %in% 
          c(VC$m3)) 
        EQ4P <- EQ5P <- EQ6P <- EQ7P <- EQ8P <- EQ9P <- 0
      if (VC$margins[2] %in% c(VC$m2) && VC$margins[3] %in% 
          c(VC$m3)) 
        EQ4P <- EQ5P <- EQ6P <- EQ7P <- EQ8P <- 0
      if (VC$margins[2] %in% c(VC$m3) && VC$margins[3] %in% 
          c(VC$m2)) 
        EQ4P <- EQ5P <- EQ6P <- EQ7P <- EQ8P <- 0
    }
    S.h <- adiag(EQ1P, EQ2P, EQ3P, EQ4P, EQ5P, EQ6P, EQ7P, 
                 EQ8P, EQ9P)
  }
  if ("kappad" %in% names(qu.mag)){
    S.h <- cbind(S.h, 0)
    S.h <- cbind(S.h, 0)
    S.h <- rbind(S.h, 0)
    S.h <- rbind(S.h, 0)
  }
  list(S.h = S.h, qu.mag = qu.mag)
}







# Get the problematic version of the function out of the package namespace
tmpfun <- get("gjrm", 
              envir = asNamespace("GJRM"))

# Make sure the new function has the environment of the old 
# function (there are possibly easier ways to do this -- I like 
# to get the old function out of the namespace to be sure I can do 
# it and am accessing what I want to access)
environment(gjrm) <- environment(tmpfun)
# Replace the old version of the function in the package namespace
# with your new version
assignInNamespace("gjrm", 
                  gjrm, ns = "GJRM")




# Get the problematic version of the function out of the package namespace
tmpfun <- get("pen", 
              envir = asNamespace("GJRM"))

# Make sure the new function has the environment of the old 
# function (there are possibly easier ways to do this -- I like 
# to get the old function out of the namespace to be sure I can do 
# it and am accessing what I want to access)
environment(pen) <- environment(tmpfun)
# Replace the old version of the function in the package namespace
# with your new version
assignInNamespace("pen", 
                  pen, ns = "GJRM")











