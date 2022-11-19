buildmer.fit <- function (p) {
	# Formula
	if (is.list(p$formula)) {
		if (is.data.frame(p$formula)) {
			p$tab <- p$formula
			if (is.null(p$dep)) stop("The 'formula' argument was specified using a buildmer terms list, but no dependent variable was specified using the 'dep' argument; please add a 'dep' argument to your buildmer() or related function call")
			p$formula <- build.formula(p$dep,p$tab,p$env)
		} else {
			stop("The 'formula' argument appears to be a list, but it does not seem to be a buildmer terms list (because those should be dataframes, which your formula isn't). The buildmer functions only work with regular formulas or with buildmer terms lists obtained from tabulate.formula(). If you got here trying to fit a multi-formula GAM, use buildcustom() to provide your own wrapper function around it - buildmer doesn't know how to manipulate mgcv's list formulas natively.")
		}
	} else {
		p$dep <- as.character(p$formula[2])
		p$tab <- tabulate.formula(p$formula)
	}

	# Include
	if (!is.null(p$include)) {
		if (inherits(p$include,'character')) {
			if ('1' %in% p$include) {
				intercept <- TRUE
				if (length(p$include) > 1) {
					# reformulate needs at least one term; it will do the right thing if this is only the intercept, but not if it's the intercept plus something else
					p$include <- p$include[p$include != '1']
				}
			} else {
				intercept <- FALSE
			}
			p$include <- stats::reformulate(p$include,intercept=intercept)
		}
		if (inherits(p$include,'formula')) {
			p$include <- tabulate.formula(p$include)
		}
	}

	# REML
	if (isTRUE(p$REML)) {
		# Force on
		p$force.reml <- TRUE
	} else if (isFALSE(p$REML)) {
		# Force off
		p$can.use.reml <- FALSE
	} else {
		# Default case, in which case one optimization can be applied:
		if (all(p$crit.name %in% c('deviance','devexp','F'))) {
			p$can.use.reml <- FALSE
			p$force.reml <- TRUE
		}
	}

	# Parallel
	if (is.null(p$cl)) {
		p$parallel <- FALSE
		p$parply <- lapply
		cleanup.cluster <- FALSE
	} else {
		p$parallel <- TRUE
		p$parply <- function (x,fun,...) parallel::parLapply(p$cl,x,fun,...)
		if (is.numeric(p$cl)) {
			p$cl <- parallel::makeCluster(p$cl,outfile='')
			cleanup.cluster <- TRUE
		} else {
			cleanup.cluster <- FALSE
		}
	}

	# Let's go
	p$reml <- p$can.use.reml || p$force.reml
	p$ordered <- ''
	crits <- p$crit
	if (length(crits) == 1) crits <- sapply(1:length(p$direction),function (i) crits)
	if (length(p$direction)) {
		if (length(crits) != length(p$direction)) {
			stop("Arguments for 'crit' and 'direction' don't make sense together -- they should have the same lengths!")
		}
		if (length(p$direction)) {
			for (i in 1:length(p$direction)) {
				p <- do.call(p$direction[[i]],list(p=within.list(p,{ crit <- crits[[i]] })))
			}
		}
		if (any(i <- names(p$results) %in% c('LRT','F'))) {
			p$results[,i] <- exp(p$results[,i])
		}
	}
	if (is.null(p$model)) {
		progress(p,'Fitting the final model')
		p$model <- p$parply(list(p),p$fit,p$formula)[[1]]
	}
	if (cleanup.cluster) {
		parallel::stopCluster(p$cl)
	}
	p
}

buildmer.finalize <- function (p) {
	ret <- mkBuildmer(model=p$model,p=p)
	ret@p$in.buildmer <- TRUE
	if (p$calc.anova) {
		ret@anova <- anova.buildmer(ret,ddf=p$ddf)
	}
	if (p$calc.summary) {
		ret@summary <- summary.buildmer(ret,ddf=p$ddf)
	}
	ret@p$in.buildmer <- FALSE
	ret
}

calcWald <- function (table,col.ef,col.df=0) {
	ef <- table[,col.ef]
	if (col.df) {
		df <- table[,col.df]
		p <- matrix(stats::pchisq(df*ef,df,lower.tail=FALSE))
		colnames(p) <- 'Pr(>F)'
	} else {
		p <- matrix(stats::pnorm(abs(ef),lower.tail=FALSE)*2)
		colnames(p) <- 'Pr(>|t|)'
	}
	cbind(table,p)
}

check.ddf <- function (model,ddf) {
	if (is.null(ddf)) {
		return('Wald')
	}
	valid <- c('Wald','lme4','Satterthwaite','Kenward-Roger','KR')
	i <- pmatch(ddf,valid)
	if (is.na(i)) {
		warning("Invalid ddf specification, possible options are 'Wald', 'lme4', 'Satterthwaite', 'Kenward-Roger'")
		ddf <- 'lme4'
	} else {
		ddf <- valid[i]
		if ('ddf' == 'KR') {
			ddf <- 'Kenward-Roger'
		}
	}

	if (ddf %in% c('Wald','lme4')) {
		return(ddf)
	} else {
		fam <- family(model)$family
		if (fam %in% c('binomial','poisson')) {
			warning('Denominator degrees of freedom do not apply to binomial or Poisson models, as those models have a known scale parameter; returning exact ddf instead')
			return('lme4')
		}
		if (startsWith(fam,'Negative Binomial')) {
			warning('Satterthwaite/Kenward-Roger denominator degrees of freedom are not available for negative binomial models; returning Wald ddf instead')
			return('lme4') #not Wald -> glmer.nb already does Wald itself
		}
		if (fam != 'gaussian' && ddf == 'Kenward-Roger') {
			warning('Kenward-Roger denominator degrees of freedom are only available for *linear* mixed models; returning Satterthwaite ddf instead')
			ddf <- 'Satterthwaite'
		}
	}

	if (!requireNamespace('lmerTest',quietly=TRUE)) {
		warning('lmerTest package is not available, could not calculate ',ddf,' denominator degrees of freedom')
		return('lme4')
	}
	if (ddf == 'Kenward-Roger' && !requireNamespace('pbkrtest',quietly=TRUE)) {
		warning('pbkrtest package is not available, could not calculate Kenward-Roger denominator degrees of freedom')
		return('lme4')
	}
	ddf
}

decompose.random.terms <- function (terms) {
	terms <- lapply(terms,function (x) {
		x <- unwrap.terms(x,inner=TRUE)
		g <- unwrap.terms(x[3])
		indep <- x[[1]] == '||'
		terms <- as.character(x[2])
		terms <- unwrap.terms(terms,intercept=TRUE)
		termlist <- if (indep) lapply(terms,identity) else list(terms)
		# We may have multiple grouping factors, e.g. in (x|a/b) constructions. We need to duplicate the termlist for each.
		ret <- rep(termlist,length(g))
		names(ret) <- rep(g,each=length(termlist))
		ret
	})
	unlist(terms,recursive=FALSE)
}

get.random.list <- function (formula) {
	bars <- lme4::findbars(formula)
	groups <- unique(sapply(bars,function (x) x[[3]]))
	randoms <- lapply(groups,function (g) {
		terms <- bars[sapply(bars,function (x) x[[3]] == g)]
		terms <- lapply(terms,function (x) x[[2]])
		terms <- lapply(terms,function (x) unravel(x,'+'))
		terms <- unique(sapply(terms,as.character))
		unique(unlist(terms))
	})
	names(randoms) <- groups
	randoms
}

has.smooth.terms <- function (formula) length(mgcv::interpret.gam(formula)$smooth.spec) > 0
is.smooth.term <- function (term) has.smooth.terms(mkForm(list(term)))
is.random.term <- function (term) {
	is.bar <- function (x) x == '|' || x == '||'
	term <- mkTerm(term)
	if (is.name(term)) return(FALSE)
	if (is.bar(term[[1]])) return(TRUE)
	if (term[[1]] == '(' && is.bar(term[[2]][[1]])) return(TRUE)
	FALSE
}
mkForm <- function (term) stats::as.formula(paste0('~',term))
mkTerm <- function (term) mkForm(term)[[2]]

progress <- function (p,...) {
	text <- sapply(list(...),function (x) as.character(list(x)))
	text <- paste0(text,collapse='')
	text <- strwrap(text,exdent=4)
	text <- paste0(text,collapse='\n')
	if (!p$quiet) {
		message(text)
	}
	text
}

unpack.smooth.terms <- function (x) {
	fm <- stats::as.formula(paste0('~',list(x)))
	if (!has.smooth.terms(fm)) return(as.character(list(x)))
	smooth.args <- fm[[2]][2:length(fm[[2]])]
	if (!all(is.null(names(smooth.args)))) smooth.args <- smooth.args[names(smooth.args) %in% c('','by')]
	unlist(lapply(smooth.args,function (x) as.character(unravel(x))))
}

unravel <- function (x,sym=c(':','interaction')) {
	if (length(x) == 1) return(as.character(x))
	if (as.character(x[[1]]) %in% sym) return(c(unravel(x[[2]],sym=sym),as.character(list(x[[3]]))))
	if (length(x) == 2) return(as.character(list(x))) #e.g.: 'scale(x)','I(365*Days)'
	# we've gotten as deep as we can go: what we now have is, e.g., :(a,:(b,c)) when sym='+'
	as.character(list(x))
}

unwrap.terms <- function (terms,inner=FALSE,intercept=FALSE) {
	form <- stats::as.formula(paste0('~',terms))
	terms <- terms(form,keep.order=TRUE)
	if (intercept) intercept <- attr(terms,'intercept')
	if (inner) return(terms[[2]])
	terms <- attr(terms,'term.labels')
	if (intercept) terms <- c('1',terms)
	terms
}
