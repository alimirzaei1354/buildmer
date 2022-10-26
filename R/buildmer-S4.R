#' The buildmer class
#' 
#' This is a simple convenience class that allows `anova' and `summary' calls to fall through to the underlying model object, while retaining buildmer's iteration history. If you need to use the final model for other things, such as prediction, access it through the `model' slot of the buildmer class object.
#' @slot model The final model containing only the terms that survived elimination
#' @slot p Parameters used during the fitting process
#' @slot anova The model's ANOVA, if the model was built with `anova=TRUE'
#' @slot summary The model's summary, if the model was built with `summary=TRUE'
#' @seealso \code{\link{buildmer}}
#' @importFrom methods new
#' @examples
#' # Manually create a bare-bones buildmer object:
#' model <- lm(Sepal.Length ~ Petal.Length,iris)
#' p <- list(in.buildmer=FALSE)
#' library(buildmer)
#' bm <- mkBuildmer(model=model,p=p,anova=NULL,summary=NULL)
#' summary(bm)
#' @export mkBuildmer
mkBuildmer <- setClass('buildmer',slots=list(model='ANY',p='list',anova='ANY',summary='ANY'))

#' @method show buildmer
#' @importFrom methods show
#' @export
show.buildmer <- function (object) {
	methods::show(object@model)
	if (length(object@p$messages)) {
		cat('\nWarning messages:\n\n')
		cat(object@p$messages)
	}
}
setMethod('show','buildmer',show.buildmer)

#' @method anova buildmer
#' @importFrom stats anova
#' @export
anova.buildmer <- function (object,...) try({
	if (length(object@p$messages)) warning(object@p$messages)
	dots <- list(...)
	ddf <- dots$ddf
	type <- dots$type

	if (!is.null(ddf) && !inherits(object@model,'merMod') && !object@p$in.buildmer) warning("Ignoring 'ddf' specification as this is not an lme4 linear mixed model")
	if (!is.null(object@anova) && is.null(ddf)) return(object@anova)
	if (inherits(object@model,'lme')) {
		if (is.null(type) || type == 3) type <- 'marginal'
		if                  (type == 1) type <- 'sequential'
		return(stats::anova(object@model,type=type))
	}
	if (!is.null(type)) {
		if (inherits(object@model,'lmerMod')) {
			if (!type %in% c(1,3)) {
				warning("Invalid 'type' argument, allowed options are 1 and 3. Resetting to type 3")
				type <- 3
			}
		}
		else warning("Ignoring 'type' argument as this is not a linear mixed model")
	}
	if (inherits(object@model,'glmmTMB')) stop('ANOVA is not available for glmmTMB fits')
	if (inherits(object@model,'MixMod')) stop('buildmer ANOVA is not available for GLMMadaptive fits; please set up an L matrix and use anova(x@model) by hand')
	if (any(names(object@model) == 'gam')) return(stats::anova(object@model$gam))
	if (!inherits(object@model,'merMod')) return(stats::anova(object@model))

	ddf <- check.ddf(object@model,ddf)
	if (is.null(type)) type <- 3
	if (ddf %in% c('Wald','lme4')) {
		table <- if (inherits(object@model,'lmerModLmerTest')) stats::anova(object@model,ddf='lme4',type=type) else stats::anova(object@model)
		if (ddf == 'Wald') {
			table <- calcWald(table,4,col.df=1)
			attr(table,'heading') <- paste('ANOVA based on type',utils::as.roman(type),'SS\n(p-values based on the Wald chi-square approximation)')
		}
	} else {
		if (!inherits(object@model,'lmerModLmerTest')) {
			object@model <- lmerTest::as_lmerModLmerTest(object@model)
		}
		table <- stats::anova(object@model,ddf=ddf,type=type)
	}
	return(table)
})

#' @method summary buildmer
#' @export
summary.buildmer <- function (object,...) try({
	if (length(object@p$messages)) warning(object@p$messages)
	dots <- list(...)
	ddf <- dots$ddf
	if (!is.null(ddf) && !inherits(object@model,'merMod') && !object@p$in.buildmer) warning("Ignoring 'ddf' specification as this is not an lme4 linear mixed model")
	if (!is.null(object@summary) && is.null(ddf)) return(object@summary)
	if (any(names(object@model) == 'gam')) return(summary(object@model$gam))
	if (!inherits(object@model,'merMod')) return(summary(object@model))

	ddf <- check.ddf(object@model,ddf)
	if (ddf %in% c('Wald','lme4')) {
		table <- if (inherits(object@model,'lmerModLmerTest')) summary(object@model,ddf='lme4') else summary(object@model)
		if (ddf == 'Wald') {
			table$coefficients <- calcWald(table$coefficients,3)
			table$methTitle <- paste0(table$methTitle,'\n(p-values based on Wald z-scores)')
		}
	} else {
		if (!inherits(object@model,'lmerModLmerTest')) object@model <- lmerTest::as_lmerModLmerTest(object@model)
		table <- summary(object@model,ddf=ddf)
	}
	return(table)
})
setMethod('summary','buildmer',summary.buildmer)

setGeneric('diag')
#' Diagonalize the random-effect covariance structure, possibly assisting convergence
#' @param x A model formula.
#' @return The formula with all random-effect correlations forced to zero, per Pinheiro & Bates (2000)
#' @examples
#' # 1. Create explicit columns for factor variables
#' library(buildmer)
#' vowels <- cbind(vowels,model.matrix(~vowel,vowels))
#' # 2. Create formula with diagonal covariance structure
#' form <- diag(f1 ~ (vowel1+vowel2+vowel3+vowel4)*timepoint*following + 
#' 	     ((vowel1+vowel2+vowel3+vowel4)*timepoint*following | participant) +
#' 	     (timepoint | word))
#' # 3. Convert formula to buildmer terms list, grouping terms starting with 'vowel'
#' terms <- tabulate.formula(form,group='vowel[^:]')
#' # 4. Directly pass the terms object to buildmer, using the 'dep' argument to specify the
#' # dependent variable
#' \donttest{model <- buildmer(terms,data=vowels,buildmerControl=list(dep='f1'))}
#' @export
setMethod('diag','formula',function (x) {
	dep <- if (length(x) < 3) '' else as.character(x[2])
	tab <- tabulate.formula(x)
	ok <- !is.na(tab$index)
	tab$index[ok] <- 1:sum(ok)
	build.formula(dep,tab,parent.frame())
})

#sapply(c('MixMod','bam','clm','clmm','gam','glm','lm','glmmTMB','gls','JuliaCall','lme','nlme','lmerMod','glmerMod','lmerModLmerTest','lmertree','glmertree','lmtree','glmtree','multinom','nnet'),function (x) methods(class=x)) %>% unlist %>% sapply(. %>% strsplit('.',fixed=T) %>% .[[1]] %>% .[1:(length(.)-1)] %>% paste0(collapse='.')) %>% unique %>% .[!endsWith(.,'-method')] %>% .[!. %in% c('anova','summary','show')] %>% sapply(function (x) {
#	forms <- names(formals(x))
#	forms2 <- paste0(forms[-1],collapse=',')
#	formsfull <- paste0(forms,collapse=',')
#	cat("#' @method",x,'buildmer\n')
#	cat("#' @export\n")
#	cat(paste0(x,'.buildmer <- function(',formsfull,') ',x,'(',forms[1],'=',forms[1],'@model,',forms2,')\n'))
#	x
#}) -> x
#' @method coef buildmer
#' @importFrom stats coef
#' @export
coef.buildmer <- function (object,...) coef(object=object@model,...)
#' @method confint buildmer
#' @importFrom stats confint
#' @export
confint.buildmer <- function (object,...) confint(object=object@model,...)
#' @method family buildmer
#' @importFrom stats family
#' @export
family.buildmer <- function (object,...) family(object=object@model,...)
#' @method fitted buildmer
#' @importFrom stats fitted
#' @export
fitted.buildmer <- function (object,...) fitted(object=object@model,...)
#' @method fixef buildmer
#' @importFrom nlme fixef
#' @export
fixef.buildmer <- function (object,...) fixef(object=object@model,...)
#' @method formula buildmer
#' @importFrom stats formula
#' @export
formula.buildmer <- function (x,...) formula(x=x@model,...)
#' @method logLik buildmer
#' @importFrom stats logLik
#' @export
logLik.buildmer <- function (object,...) logLik(object=object@model,...)
#' @method model.frame buildmer
#' @importFrom stats model.frame
#' @export
model.frame.buildmer <- function (formula,...) model.frame(formula=formula@model,...)
#' @method model.matrix buildmer
#' @importFrom stats model.matrix
#' @export
model.matrix.buildmer <- function (object,...) model.matrix(object=object@model,...)
#' @method nobs buildmer
#' @importFrom stats nobs
#' @export
nobs.buildmer <- function (object,...) nobs(object=object@model,...)
#' @method predict buildmer
#' @importFrom stats predict
#' @export
predict.buildmer <- function (object,...) predict(object=object@model,...)
#' @method print buildmer
#' @export
print.buildmer <- function (x,...) print(x=x@model,...)
#' @method ranef buildmer
#' @importFrom nlme ranef
#' @export
ranef.buildmer <- function (object,...) ranef(object=object@model,...)
#' @method residuals buildmer
#' @importFrom stats residuals
#' @export
residuals.buildmer <- function (object,...) residuals(object=object@model,...)
#' @method simulate buildmer
#' @importFrom stats simulate
#' @export
simulate.buildmer <- function (object,...) simulate(object=object@model,...)
#' @method terms buildmer
#' @importFrom stats terms
#' @export
terms.buildmer <- function (x,...) terms(x=x@model,...)
#' @method vcov buildmer
#' @importFrom stats vcov
#' @export
vcov.buildmer <- function (object,...) vcov(object=object@model,...)
#' @method cooks.distance buildmer
#' @importFrom stats cooks.distance
#' @export
cooks.distance.buildmer <- function (model,...) cooks.distance(model=model@model,...)
#' @method influence buildmer
#' @importFrom stats influence
#' @export
influence.buildmer <- function (model,...) influence(model=model@model,...)
#' @method plot buildmer
#' @importFrom graphics plot
#' @export
plot.buildmer <- function (x,...) plot(x=x@model,...)
#' @method add1 buildmer
#' @importFrom stats add1
#' @export
add1.buildmer <- function (object,...) add1(object=object@model,...)
#' @method deviance buildmer
#' @importFrom stats deviance
#' @export
deviance.buildmer <- function (object,...) deviance(object=object@model,...)
#' @method drop1 buildmer
#' @importFrom stats drop1
#' @export
drop1.buildmer <- function (object,...) drop1(object=object@model,...)
#' @method effects buildmer
#' @importFrom stats effects
#' @export
effects.buildmer <- function (object,...) effects(object=object@model,...)
#' @method extractAIC buildmer
#' @importFrom stats extractAIC
#' @export
extractAIC.buildmer <- function (fit,...) extractAIC(fit=fit@model,...)
#' @method profile buildmer
#' @importFrom stats profile
#' @export
profile.buildmer <- function (fitted,...) profile(fitted=fitted@model,...)
#' @method rstandard buildmer
#' @importFrom stats rstandard
#' @export
rstandard.buildmer <- function (model,...) rstandard(model=model@model,...)
#' @method rstudent buildmer
#' @importFrom stats rstudent
#' @export
rstudent.buildmer <- function (model,...) rstudent(model=model@model,...)
#' @method weights buildmer
#' @importFrom stats weights
#' @export
weights.buildmer <- function (object,...) weights(object=object@model,...)
#' @method alias buildmer
#' @importFrom stats alias
#' @export
alias.buildmer <- function (object,...) alias(object=object@model,...)
#' @method case.names buildmer
#' @importFrom stats case.names
#' @export
case.names.buildmer <- function (object,...) case.names(object=object@model,...)
#' @method dfbeta buildmer
#' @importFrom stats dfbeta
#' @export
dfbeta.buildmer <- function (model,...) dfbeta(model=model@model,...)
#' @method dfbetas buildmer
#' @importFrom stats dfbetas
#' @export
dfbetas.buildmer <- function (model,...) dfbetas(model=model@model,...)
#' @method dummy.coef buildmer
#' @importFrom stats dummy.coef
#' @export
dummy.coef.buildmer <- function (object,...) dummy.coef(object=object@model,...)
#' @method hatvalues buildmer
#' @importFrom stats hatvalues
#' @export
hatvalues.buildmer <- function (model,...) hatvalues(model=model@model,...)
#' @method kappa buildmer
#' @export
kappa.buildmer <- function (z,...) kappa(z=z@model,...)
#' @method labels buildmer
#' @export
labels.buildmer <- function (object,...) labels(object=object@model,...)
#' @method proj buildmer
#' @importFrom stats proj
#' @export
proj.buildmer <- function (object,...) proj(object=object@model,...)
#' @method qqnorm buildmer
#' @importFrom stats qqnorm
#' @export
qqnorm.buildmer <- function (y,...) qqnorm(y=y@model,...)
#' @method qr buildmer
#' @export
qr.buildmer <- function (x,...) qr(x=x@model,...)
#' @method variable.names buildmer
#' @importFrom stats variable.names
#' @export
variable.names.buildmer <- function (object,...) variable.names(object=object@model,...)
#' @method df.residual buildmer
#' @importFrom stats df.residual
#' @export
df.residual.buildmer <- function (object,...) df.residual(object=object@model,...)
#' @method sigma buildmer
#' @importFrom stats sigma
#' @export
sigma.buildmer <- function (object,...) sigma(object=object@model,...)
#' @method VarCorr buildmer
#' @importFrom nlme VarCorr
#' @export
VarCorr.buildmer <- function (x,...) VarCorr(x=x@model,...)
#' @method ACF buildmer
#' @importFrom nlme ACF
#' @export
ACF.buildmer <- function (object,...) ACF(object=object@model,...)
#' @method augPred buildmer
#' @importFrom nlme augPred
#' @export
augPred.buildmer <- function (object,...) augPred(object=object@model,...)
#' @method comparePred buildmer
#' @importFrom nlme comparePred
#' @export
comparePred.buildmer <- function (object1,...) comparePred(object1=object1@model,...)
#' @method getData buildmer
#' @importFrom nlme getData
#' @export
getData.buildmer <- function (object) getData(object=object@model)
#' @method getGroups buildmer
#' @importFrom nlme getGroups
#' @export
getGroups.buildmer <- function (object,form,level,data,sep) getGroups(object=object@model,form,level,data,sep)
#' @method getGroupsFormula buildmer
#' @importFrom nlme getGroupsFormula
#' @export
getGroupsFormula.buildmer <- function (object,asList,sep) getGroupsFormula(object=object@model,asList,sep)
#' @method getResponse buildmer
#' @importFrom nlme getResponse
#' @export
getResponse.buildmer <- function (object,form) getResponse(object=object@model,form)
#' @method getVarCov buildmer
#' @importFrom nlme getVarCov
#' @export
getVarCov.buildmer <- function (obj,...) getVarCov(obj=obj@model,...)
#' @method intervals buildmer
#' @importFrom nlme intervals
#' @export
intervals.buildmer <- function (object,...) intervals(object=object@model,...)
#' @method update buildmer
#' @importFrom stats update
#' @export
update.buildmer <- function (object,...) update(object=object@model,...)
#' @method Variogram buildmer
#' @importFrom nlme Variogram
#' @export
Variogram.buildmer <- function (object,...) Variogram(object=object@model,...)
#' @method pairs buildmer
#' @importFrom graphics pairs
#' @export
pairs.buildmer <- function (x,...) pairs(x=x@model,...)
#' @method step buildmer
#' @importFrom stats step
#' @export
step.buildmer <- function (object,...) step(object=object@model,...)
#' @method getME buildmer
#' @importFrom lme4 getME
#' @export
getME.buildmer <- function (object,...) getME(object=object@model,...)
#' @method isLMM buildmer
#' @importFrom lme4 isLMM
#' @export
isLMM.buildmer <- function (x,...) isLMM(x=x@model,...)
#' @method refit buildmer
#' @importFrom lme4 refit
#' @export
refit.buildmer <- function (object,...) refit(object=object@model,...)
