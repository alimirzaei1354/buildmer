callfixup <- function (p,substitute.fun,call,patch.family) {
	for (x in intersect(NSENAMES,names(p$args))) {
		call[x] <- p$call$args[x]
	}
	call[[1]] <- substitute.fun
	call$data <- p$call$data
	if (patch.family) {
		call$family <- p$call$family
	}
	call
}
run <- function (fun,args,quiet) {
	if (quiet) {
		suppressMessages(suppressWarnings(try(do.call(fun,args),silent=TRUE)))
	} else {
		suppressWarnings(try(do.call(fun,args)))
	}
}

patch.GLMMadaptive <- function (p,fun,args) {
	model <- run(fun,args,p$quiet)
	if (inherits(model,'try-error')) {
		return(model)
	}
	model$call <- callfixup(p,substitute(fun),model$call,TRUE)
	model
}

patch.gamm <- function (p,fun,args) {
	model <- run(fun,args,p$quiet)
	if (inherits(model,'try-error')) {
		return(model)
	}
	model$lme$call <- callfixup(p,substitute(fun),model$lme$call,TRUE)
	model
}

patch.gamm4 <- function (p,fun,args) {
	model <- run(fun,args,p$quiet)
	if (inherits(model,'try-error')) {
		return(model)
	}
	model$mer@call <- callfixup(p,substitute(fun),model$mer@call,TRUE)
	model
}

patch.lm <- function (p,fun,args) {
	model <- run(fun,args,p$quiet)
	if (inherits(model,'try-error')) {
		return(model)
	}
	model$call <- callfixup(p,substitute(fun),model$call,!p$is.gaussian)
	model
}

patch.lmer <- function (p,fun,args) {
	model <- run(fun,args,p$quiet)
	if (inherits(model,'try-error')) {
		return(model)
	}
	model@call <- callfixup(p,substitute(fun),model@call,!p$is.gaussian)
	model
}

patch.mertree <- function (p,fun,args) {
	model <- run(fun,args,p$quiet)
	if (inherits(model,'try-error')) {
		return(model)
	}
	eltname <- if (p$is.gaussian) 'lmer' else 'glmer'
	if (!converged(model[[eltname]],p$singular.ok,p$grad.tol,p$hess.tol)) {
		return(model[[eltname]])
	}
	model$call <- callfixup(p,substitute(fun),model$call,!p$is.gaussian)
	model[[eltname]]@call <- callfixup(p,substitute(fun),model[[eltname]]@call,!p$is.gaussian)
	model$call$ctrl <- p$call$args$control
	model[[eltname]]@call$control <- if (p$is.gaussian) p$call$args$lmer.control else p$call$args$glmer.control
	model
}
