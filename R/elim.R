get2LL <- function (m) as.numeric(-2*stats::logLik(m))
getdevexp <- function (m) {
	if (all(c('deviance','null.deviance') %in% names(m))) return(1-m$deviance/m$null.deviance)
	ff <- fitted(m)
	rr <- stats::resid(m)
	stats::cor(ff,ff+rr)^2
}
safeNdf <- function (m) {
	if (inherits(m,c('nnet','multinom'))) {
		# no nobs method
		attr(logLik(m),'df')
	} else {
		# for some models, e.g. GAMs, this is safer than attr(logLik(m),'df')
		nobs(m) - safeRdf(m)
	}
}
safeRdf <- function (m) {
	if (inherits(m,c('nnet','multinom'))) {
		# no df.residual method nor nobs method
		nrow(fitted(m)) - attr(logLik(m),'df')
	} else {
		if (is.null(rdf <- df.residual(m))) {
			# apparently, no df.residual implementation; if there is a nobs method, we can work around this
			# (this is the case for at least MixMod, clm, clmm, gls, and lme)
			rdf <- nobs(m) - attr(logLik(m),'df')
		}
		rdf
	}
}

crit.AIC <- function (p,ref,alt) if (is.null(ref)) stats::AIC(alt) else stats::AIC(alt) - stats::AIC(ref)
crit.BIC <- function (p,ref,alt) if (is.null(ref)) stats::BIC(alt) else stats::BIC(alt) - stats::BIC(ref)
crit.F <- function (p,ref,alt) {
	r2_alt  <- getdevexp(alt)
	ddf_alt <- safeRdf(alt)
	ndf_alt <- safeNdf(alt)
	if (is.null(ref)) {
		r2_ref <- ndf_ref <- 0
	} else {
		r2_ref  <- getdevexp(ref)
		ndf_ref <- safeNdf(ref)
	}
	if (is.null(r2_alt) || is.null(r2_ref)) {
		stop('r^2 not available for this family, cannot compute the F criterion!')
	}
	Fval <- ddf_alt * (r2_alt - r2_ref) / (1 - r2_alt)
	ndf  <- ndf_alt - ndf_ref
	ddf  <- ddf_alt
	if (is.na(Fval)) {
		return(log(1))
	}
	if (Fval <= 0 || ndf <= 0) {
		return(log1p(abs(Fval))) #gives the order step some idea of which model is the least unhelpful
	}
	scale.est <- if (is.na(p$scale.est)) {
		if ('scale.estimated' %in% names(alt)) {
			alt$scale.estimated
		} else {
			TRUE #not special-cased, not binomial or poisson --> true by default
		}
	} else {
		p$scale.est
	}
	if (scale.est) {
		stats::pf(Fval,ndf,ddf,lower.tail=FALSE,log.p=TRUE)
	} else {
		stats::pchisq(ndf*Fval,ndf,lower.tail=FALSE,log.p=TRUE)
	}
}
crit.LRT <- function (p,ref,alt) {
	if (is.null(ref)) {
		chLL <- get2LL(alt)
		chdf <- safeNdf(alt)
		f1   <- ~0
	} else {
		chLL <- get2LL(ref) - get2LL(alt)
		chdf <- safeRdf(ref) - safeRdf(alt)
		f1   <- formula(ref)
	}
	if (chdf <= 0) {
		return(0)
	}

	# If the two models differ in lme4-style random effects, we need to correct the p-value
	# We cannot use the formula stored in p here because that may still go through re2mgcv
	f2 <- formula(alt)
	if (inherits(f1,'formula') && inherits(f2,'formula')) {
		# may be lme4 formulas
		tab.ref <- tabulate.formula(f1)
		tab.alt <- tabulate.formula(f2)
		if (diff.re <- sum(!is.na(tab.alt$grouping)) - sum(!is.na(tab.ref$grouping))) {
			# Stram & Lee (1994): mixture of chisq(chdf) and chisq(chdf-1)
			p1 <- stats::pchisq(chLL,chdf        ,lower.tail=FALSE,log.p=TRUE) - log(2)
			p2 <- stats::pchisq(chLL,chdf-diff.re,lower.tail=FALSE,log.p=TRUE) - log(2)
			return(log(exp(p1) + exp(p2)))
		}
	}
	stats::pchisq(chLL,chdf,lower.tail=FALSE,log.p=TRUE)
}
crit.2LL <- function (p,ref,alt) if (is.null(ref)) get2LL(alt) else get2LL(alt) - get2LL(ref)
crit.LL <- crit.2LL
crit.devexp <- function (p,ref,alt) if (is.null(ref)) getdevexp(alt) else getdevexp(alt) - getdevexp(ref)
crit.deviance <- crit.devexp

elim.AIC <- function (diff) diff > -.001
elim.BIC <- elim.AIC
elim.F   <- function (logp) exp(logp) >= .05
elim.LRT <- elim.F
elim.2LL <- elim.AIC
elim.LL  <- elim.AIC
elim.devexp <- elim.AIC
elim.deviance <- elim.AIC
