"locfit" <- 
function(formula, data = sys.frame(sys.parent()), weights = 1, cens = NULL, base = 0, 
	subset, geth = F, ..., lfproc = locfit.raw)
{
	Terms <- terms(formula, data = data, specials = c("left", "right", 
		"ang", "lp"))
	attr(Terms, "intercept") <- 0
	m <- match.call()
	m[[1]] <- as.name("model.frame")
	z <- pmatch(names(m), c("formula", "data", "weights", "cens", "base", 
		"subset"))
	for(i in length(z):2)
		if(is.na(z[i])) m[[i]] <- NULL
	frm <- eval(m, sys.frame(sys.parent()))
	vnames <- as.character(attributes(Terms)$variables)[-1]
	if(attr(Terms, "response")) {
		y <- model.extract(frm, response)
		yname <- deparse(formula[[2]])
		vnames <- vnames[-1]
	}
	else {
		y <- yname <- NULL
	}
	x <- as.matrix(frm[, vnames])
	d <- ncol(x)
	sty <- rep(1, d)
	spec <- attr(Terms, "specials")
	if(!is.null(spec)) {
		nr <- attr(Terms, "response")
		sty[spec$left - nr] <- 5
		sty[spec$right - nr] <- 6
		sty[spec$ang - nr] <- 3
	}
	if(!missing(weights))
		weights <- model.extract(frm, weights)
	if(!missing(cens))
		cens <- model.extract(frm, cens)
	if(!missing(base))
		base <- model.extract(frm, base)
	ret <- lfproc(x, y, weights = weights, cens = cens, base = base, geth
		 = geth, sty = sty, ...)
	if(any(geth == 1:4))
		return(ret)
	ret$terms <- Terms
	ret$call <- match.call()
	if(!is.null(yname))
		ret$yname <- yname
	ret$vnames <- vnames
	ret$frame <- sys.frame(sys.parent())
	ret
}
"locfit.raw" <- 
function(x, y, weights = 1, cens = NULL, base = 0, subset, xlim, flim, scale = 
	rep(1, d), alpha = 0.7, ev = "tree", deg = 2, family, link = "default", 
	maxk = 50, kern = "tcub", kt = "sph", itype = "default", acri = "cp", 
	mint = 20, maxit = 20, cut = 0.8, dc = F, geth = F, renorm = F, mg = 10,
	deriv = numeric(0), sty = rep(1, d))
{
	if(!is.matrix(x)) {
		vnames <- deparse(substitute(x))
		x <- matrix(x, ncol = 1)
		d <- 1
	}
	else {
		d <- ncol(x)
		if(is.null(dimnames(x)))
			vnames <- paste("x", 1:d, sep = "")
		else vnames <- dimnames(x)[[2]]
	}
	n <- nrow(x)
	if((!missing(y)) && (!is.null(y))) {
		yname <- deparse(substitute(y))
		if(missing(family))
			family <- if(is.logical(y)) "binomial" else "qgaussian"
	}
	else {
		if(missing(family))
			family <- if(is.null(cens)) "density" else "hazard"
		y <- 0
		yname <- family
	}
	if(is.null(cens))
		cens <- 0
	kt <- pmatch(kt, c("sph", "prod", "ang", "center", "left", "right"))
	if(any(kt == c(3, 5, 6)))
		stop("kt=ang,left,right no longer used. Use eg. y~ang(x) formula"
			)
	if(length(deg) >= 2) {
		deg0 <- deg[1]
		deg <- deg[2]
	}
	else deg0 <- deg
	p <- ifelse(kt == 2, d * deg + 1, ifelse(deg == 0, 1, prod(d + (1:deg))/ prod(1:deg)))
	if(is.numeric(ev)) {
		xev <- ev
		vc <- ncm <- 0
		nvm <- length(xev)/d
		ev <- 8
	}
	else {
		ev <- pmatch(ev, c("tree", "phull", "data", "grid", "kdtree", 
			"kdcenter", "crossval"))
		switch(ev,
			{
				vc <- 2^d
				nvm <- ((maxk + 2) * vc)/2
				ncm <- 1 + 2 * maxk
			}
			,
			{
				vc <- d + 1
				nvm <- ncm <- maxk * d
			}
			,
			{
				vc <- ncm <- 0
				nvm <- n
			}
			,
			{
				vc <- 2^d
				mg <- rep(mg, length = d)
				nvm <- max(prod(mg), maxk)
				ncm <- 0
			}
			,
			{
				vc <- 2^d
				fc <- floor((cut * n * min(alpha, 1))/4)
				k <- floor((2 * n)/fc)
				ncm <- 2 * k + 1
				nvm <- ((k + 2) * vc)/2
			}
			,
			{
				vc <- 1
				fc <- floor(n * alpha)
				nvm <- 1 + floor((2 * n)/fc)
				ncm <- 2 * nvm + 1
			}
			,
			{
				vc <- ncm <- 0
				nvm <- n
			}
			)
	}
	xl <- rep(0, 2 * d)
	lset <- 0
	if(!missing(xlim)) {
		xl <- lflim(xlim, vnames, xl)
		lset <- 1
	}
	fl <- rep(0, 2 * d)
	if(!missing(flim))
		fl <- lflim(flim, vnames, fl)
	mi <- c(n, p, deg0, deg, d, 0, 0, kt, 0, mint, maxit, renorm, ev, 0, 0, 
		dc, maxk, 1, geth)
	if(any(is.na(mi)))
		print(mi)
	alpha <- c(alpha, 0, 0, 0)[1:3]
	if(length(scale) == 1)
		scale <- rep(scale, d)
	dp <- c(alpha, cut, 0, 0, 0, 0)
	w1 <- n * (p + 5) + 4 * p * p + 6 * p
	nnl <- p + d
	if(ev == 8)
		w2 <- c(xev, numeric(nvm * (d + p + 6 + nnl) + ncm))
	else w2 <- numeric(nvm * (2 * d + p + 6 + nnl) + ncm)
	w3 <- n + ncm * vc + 3 * max(nvm, ncm)
	z <- .C("Sfit",
		x = as.numeric(x),
		y = as.numeric(rep(y, length.out = n)),
		cens = as.numeric(rep(cens, length.out = n)),
		w = as.numeric(rep(weights, length.out = n)),
		base = as.numeric(rep(base, length.out = n)),
		lim = as.numeric(c(xl, fl)),
		mi = as.integer(mi),
		dp = as.numeric(dp),
		strings = c(kern, family, link, itype, acri),
		scale = as.numeric(scale),
		wdes = numeric(w1),
		wtre = w2,
		nvc = as.integer(c(nvm, ncm, nnl, 0, 0)),
		iwork = integer(w3),
		lw = as.integer(c(w1, length(w2), w3)),
		mg = as.integer(mg),
		L = numeric(c(1, n * nvm, 2 * n * (d * d + d + 1), 1, 1)[geth + 
			1]),
		kap = numeric(c(3, 3, 3, length(deriv), 3)[geth + 1]),
		deriv = as.integer(deriv),
		nd = as.integer(length(deriv)),
		sty = as.integer(sty))
	nvc <- z$nvc
	names(nvc) <- c("nvm", "ncm", "nnl", "nv", "nc")
	nvm <- nvc["nvm"]
	ncm <- nvc["ncm"]
	nv <- nvc["nv"]
	nc <- nvc["nc"]
	if(geth == 1)
		return(matrix(z$L[1:(nv * n)], ncol = nv))
	if(geth == 2)
		return(list(const = z$kap, d = d))
	if(geth == 3)
		return(z$kap)
	dp <- z$dp
	mi <- z$mi
	names(mi) <- c("n", "p", "deg0", "deg", "d", "acri", "ker", "kt", "it", 
		"mint", "mxit", "renorm", "ev", "tg", "link", "dc", "mk", "wh", 
		"geth")
	names(dp) <- c("nnalph", "fixh", "adpen", "cut", "lk", "t0", "t1", "rv"
		)
	if(geth == 4)
		return(list(residuals = z$y, var = z$wdes[n * (p + 2) + p * p + (
			1:n)], nl.df = dp["t0"] - 2))
	trans <- switch(mi["link"] - 2,
		function(x)
		x,
		exp,
		expit,
		function(x)
		1/x,
		function(x)
		x * x)
	t1 <- z$wtre
	t2 <- z$iwork[ - (1:n)]
	eva <- list(xev = t1[1:(d * nv)], coef = matrix(t1[d * nvm + 1:((p + d + 
		nnl + 6) * nvm)], nrow = nvm)[1:nv,  ], scale = z$scale)
	if(nv == 1)
		eva$coef <- matrix(eva$coef, nrow = 1)
	class(eva) <- "lfeval"
	if(nc == 0) {
		cell <- list(sv = integer(0), ce = integer(0), s = integer(0), 
			lo = as.integer(rep(0, nv)), hi = as.integer(rep(0, nv)
			))
	}
	else {
		mvc <- max(nv, nc)
		mvcm <- max(nvm, ncm)
		cell <- list(sv = t1[nvm * (2 * d + p + 6 + nnl) + 1:nc], ce = 
			t2[1:(vc * nc)], s = t2[vc * ncm + 1:mvc], lo = t2[vc * 
			ncm + mvcm + 1:mvc], hi = t2[vc * ncm + 2 * mvcm + 1:
			mvc])
	}
	ret <- list(eva = eva, cell = cell, terms = NULL, nvc = nvc, box = z$
		lim[2 * d + 1:(2 * d)], sty = sty, mi = mi, mg = z$mg, dp = dp, 
		trans = trans, vnames = vnames, yname = yname, call = 
		match.call(), frame = sys.frame(sys.parent()))
	class(ret) <- "locfit"
	ret
}
"gam.lf" <- 
function(x, y, w, xeval, ...)
{
	if(!missing(xeval)) {
		fit <- locfit.raw(x, y, weights = w, ...)
		return(predict(fit, x))
	}
	locfit.raw(x, y, weights = w, geth = 4, ...)
}
"lf" <- 
function(..., alpha = 0.7, deg = 2, scale = 1, ev = "tree", maxk = 100)
{
	if(!any(gam.slist == "lf"))
		warning("gam.slist does not include \"lf\" -- fit will be incorrect"
			)
	scall <- deparse(sys.call())
	vars <- list(...)
	x <- vars[[1]]
	if(length(vars) > 1)
		for(i in 2:length(vars))
			x <- cbind(x, vars[[i]])
	d <- if(is.matrix(x)) ncol(x) else 1
	attr(x, "alpha") <- alpha
	attr(x, "deg") <- deg
	attr(x, "scale") <- scale
	attr(x, "ev") <- ev
	attr(x, "maxk") <- maxk
	attr(x, "ncols") <- d
	real.call <- substitute(gam.lf(data[[scall]], z, w, alpha = alpha, deg
		 = deg, scale = scale, maxk = maxk, ev = ev))
	attr(x, "call") <- real.call
	attr(x, "class") <- "smooth"
	x
}
"gam.slist" <- 
c("s", "lo", "random", "lf")
"left" <- 
function(x)
x
"right" <- 
function(x)
x
"ang" <- 
function(x)
x
"fitted.locfit" <- 
function(object, data = NULL, what = "coef", se.fit = F, cv = F, tr)
{
	if(missing(data)) {
		data <- if(is.null(object$call$data)) sys.frame(sys.parent()) else eval(
				object$call$data)
	}
	if(missing(tr))
		tr <- if(what == "coef") object$trans else function(x)
			x
	mm <- locfit.matrix(object, data = data)$x
	ev <- object$mi["ev"]
	if(any(ev == c(3, 7))) {
		fit <- knots(object, what = what, tr = function(x)
		x)
		if(se.fit)
			se <- knots(object, what = "se")
	}
	else {
		pred <- preplot(object, mm, se.fit = se.fit, what = what)
		fit <- pred$fit
		if(se.fit)
			se <- pred$se.fit
	}
	if((cv) && (ev != 7)) {
		infl <- preplot(object, mm, what = "infl")$fit * residuals(
			object, type = "ldot")
		fit <- fit - infl
	}
	fit <- tr(fit)
	if(se.fit)
		fit <- cbind(fit, se)
	fit
}
"formula.locfit" <- 
function(object)
object$call$formula
"lines.locfit" <- 
function(x, m = 100, ...)
{
	newx <- lfmarg(x, m = m)[[1]]
	y <- predict(x, newx)
	lines(newx, y, ...)
}
"plot.locfit" <- 
function(x, xlim, pv, tv, mpv, mtv = 6, deriv = numeric(0), se.fit = F, tr = 
	NULL, what = "coef", get.data = F, f3d = (d == 2) && (length(tv) > 0), 
	...)
{
	d <- x$mi["d"]
	if(missing(pv))
		pv <- if(d == 1) 1 else c(1, 2)
	if(is.character(pv))
		pv <- match(pv, x$vnames)
	if(missing(tv))
		tv <- (1:d)[ - pv]
	if(is.character(tv))
		tv <- match(tv, x$vnames)
	vrs <- c(pv, tv)
	if(any(duplicated(vrs)))
		warning("Duplicated variables in pv, tv")
	if(any((vrs <= 0) | (vrs > d)))
		stop("Invalid variable numbers in pv, tv")
	if(missing(mpv))
		mpv <- if(d == 1) 100 else 40
	m <- rep(mpv, d)
	m[tv] <- mtv
	xl <- x$box
	if(!missing(xlim))
		xl <- lflim(xlim, x$vnames, xl)
	marg <- lfmarg(xl, m)
	pred <- preplot.locfit(x, marg, se.fit = se.fit, deriv = deriv, tr = tr,
		what = what, get.data = get.data, f3d = f3d)
	plot(pred, pv = pv, tv = tv, ...)
}
"points.locfit" <- 
function(x, tr, ...)
{
	d <- x$mi["d"]
	p <- x$mi["p"]
	nv <- x$nvc["nv"]
	if(d == 1) {
		if(missing(tr))
			tr <- x$trans
		x1 <- x$eva$xev
		x2 <- x$eva$coef[1:nv]
		points(x1, tr(x2), ...)
	}
	if(d == 2) {
		xx <- knots(x, what = "x")
		points(xx[, 1], xx[, 2], ...)
	}
}
"preplot.locfit" <- 
function(object, newdata = NULL, se.fit = F, deriv = numeric(0), tr = NULL, 
	what = "coef", get.data = F, f3d = F)
{
	gr <- F
	mi <- object$mi
	dim <- mi["d"]
	if(is.null(newdata) | (object$mi["ev"] == 3))
		xev <- matrix(object$eva$xev, ncol = dim, byrow = T)
	else {
		if(is.data.frame(newdata))
			xev <- as.matrix(model.frame(delete.response(object$
				terms), newdata))
		else if(is.list(newdata)) {
			gr <- T
			xev <- as.matrix(expand.grid(newdata))
		}
		else xev <- as.matrix(newdata)
	}
	m <- nrow(xev)
	mi["wh"] <- pmatch(what, c("coef", "infl", "nlx", "band", "degr"))
	if(length(deriv) > 0) {
		if(any((deriv <= 0) | (deriv > dim)) | (length(deriv) > mi[
			"deg"]))
			stop("Invalid deriv argument")
		if((is.null(tr)) | missing(tr))
			tr <- function(x)
			x
	}
	else {
		if(is.null(tr)) {
			if(mi["wh"] == 1)
				tr <- object$trans
			else tr <- function(x)
				x
		}
	}
	nr <- 1 + se.fit
	z <- .C("Spred",
		xev = as.numeric(object$eva$xev),
		coef = as.numeric(object$eva$coef),
		sv = as.numeric(object$cell$sv),
		ce = as.integer(c(object$cell$ce, object$cell$s, object$cell$lo,
			object$cell$hi)),
		x = as.numeric(xev),
		y = numeric(m * nr),
		cse = as.integer(se.fit),
		scale = as.numeric(object$eva$scale),
		m = as.integer(m),
		nvc = as.integer(object$nvc),
		mi = as.integer(mi),
		dp = as.numeric(object$dp),
		mg = as.integer(object$mg),
		deriv = as.integer(deriv),
		nd = as.integer(length(deriv)),
		sty = as.integer(object$sty))
	fhat <- z$y[1:m]
	fhat[fhat == 0.1278433] <- NA
	if(se.fit)
		sse <- sqrt(object$dp["rv"]) * z$y[m + 1:m]
	else sse <- numeric(0)
	if(!gr)
		newdata <- list(xev = xev)
	data <- if(get.data) locfit.matrix(object) else list()
	if((f3d) | (dim > 3))
		dim <- 3
	ret <- list(xev = newdata, fit = fhat, se.fit = sse, residual.scale = 
		sqrt(object$dp["rv"]), trans = tr, vnames = object$vnames, 
		yname = object$yname, dim = as.integer(dim), data = data)
	class(ret) <- "locfitpred"
	ret
}
"predict.locfit" <- 
function(object, newdata = NULL, se.fit = F, ...)
{
	pred <- preplot.locfit(object, newdata, se.fit = se.fit, ...)
	fit <- pred$trans(pred$fit)
	if(!se.fit)
		return(fit)
	return(list(fit = fit, se.fit = pred$se.fit, residual.scale = pred$
		residual.scale))
}
"print.locfit" <- 
function(x, ...)
{
	print(knots(x, ...))
	invisible(x)
}
"residuals.locfit" <- 
function(object, data = NULL, type = "deviance", studentize = F)
{
	if(missing(data)) {
		data <- if(is.null(object$call$data)) sys.frame(sys.parent()) else eval(
				object$call$data)
	}
	frm <- locfit.matrix(object, data = data)
	x <- frm$x
	y <- frm$y
	if(any(object$mi["ev"] == c(3, 7))) {
		fv <- knots(object, what = "coef", tr = function(x)
		x)
	}
	else fv <- preplot(object, x)$fit
	n <- length(y)
	ty <- pmatch(type, c("deviance", "pearson", "response", "ldot", "mean", 
		"d2", "ldd"))
	if(is.na(ty))
		stop(paste("Invalid residual type", type))
	resid <- .C("compres",
		y = as.numeric(y),
		c = as.integer(rep(frm$ce, length.out = n)),
		w = as.numeric(rep(frm$w, length.out = n)),
		th = as.numeric(fv),
		mi = as.integer(object$mi),
		ty = as.integer(ty),
		m = as.integer(n))$th
	if(studentize) {
		pinfl <- summary(preplot(object, x, what = "infl"))
		pnlx <- summary(preplot(object, x, what = "nlx"))
		resid <- resid/sqrt(1 - 2 * frm$w * pinfl + frm$w * pnlx * pnlx
			)
	}
	resid
}
"plot.locfit.1d" <- 
function(x, pv, tv, main = "", xlab = "default", ylab = x$yname, type = "l", 
	ylim, col = 1, crit = 1.96, add = F, ...)
{
	y <- x$fit
	nos <- !is.na(y)
	xev <- x$xev[[1]][nos]
	y <- y[nos]
	ord <- order(xev)
	if(xlab == "default")
		xlab <- x$vnames
	tr <- x$trans
	yy <- tr(y)
	if(length(x$se.fit) > 0) {
		if((!is.null(class(crit))) && (class(crit) == "kappa"))
			crit <- crit$crit.val
		cup <- tr((y + crit * x$se.fit))[ord]
		clo <- tr((y - crit * x$se.fit))[ord]
	}
	if(missing(ylim)) {
		if(length(x$se.fit) > 0)
			ylim <- c(min(clo), max(cup))
		else ylim <- range(yy)
		if((length(x$data) > 0) && (!is.null(x$data$y)))
			ylim <- range(c(ylim, x$data$y/x$data$sc))
	}
	if(!add) {
		plot(xev[ord], yy[ord], type = "n", xlab = xlab, ylab = ylab, 
			main = main, xlim = range(x$xev[[1]]), ylim = ylim, ...
			)
	}
	lines(xev[ord], yy[ord], type = type, col = col, ...)
	if(length(x$se.fit) > 0) {
		lines(xev[ord], cup, lty = 2)
		lines(xev[ord], clo, lty = 2)
	}
	if(length(x$data) > 0) {
		xd <- x$data$x
		yd <- x$data$y/x$data$sc
		cd <- rep(x$data$ce, length.out = length(xd))
		if(length(yd) < 2) {
			rug(xd[cd == 0])
			if(any(cd == 1))
				rug(xd[cd == 1], ticksize = 0.015)
		}
		else {
			plotbyfactor(xd, yd, cd, col = col, pch = c("o", "+"), 
				add = T)
		}
	}
	invisible(NULL)
}
"plot.locfit.2d" <- 
function(x, pv, tv, main, type = "contour", xlab, ...)
{
	if(!is.list(x$xev))
		stop("Can only plot from grids")
	if(missing(xlab))
		xlab <- x$vnames
	tr <- x$trans
	m1 <- x$xev[[1]]
	m2 <- x$xev[[2]]
	y <- matrix(tr(x$fit))
	if(type == "contour")
		contour(m1, m2, matrix(y, nrow = length(m1)), xlab = xlab[1], 
			ylab = xlab[2], ...)
	if(type == "image")
		image(m1, m2, matrix(y, nrow = length(m1)), xlab = xlab[1], 
			ylab = xlab[2], ...)
	if((length(x$data) > 0) & any(type == c("contour", "image"))) {
		xd <- x$data$x
		ce <- rep(x$data$ce, length = nrow(xd))
		points(xd[ce == 0, 1], xd[ce == 0, 2], pch = "o")
		if(any(ce == 1))
			points(xd[ce == 1, 1], xd[ce == 1, 2], pch = "+")
	}
	if(type == "persp") {
		nos <- is.na(y)
		y[nos] <- min(y[!nos])
		persp(m1, m2, matrix(y, nrow = length(m1)), xlab = xlab[1], 
			ylab = xlab[2], ...)
	}
	if(!missing(main))
		title(main = main)
	invisible(NULL)
}
"plot.locfit.3d" <- 
function(x, main = "", pv, tv, type = "level", pred.lab = x$vnames, resp.lab = 
	x$yname, crit = 1.96, ...)
{
	if(!is.list(x$xev))
		stop("Can only plot from grids")
	newx <- as.matrix(expand.grid(x$xev))
	newy <- x$trans(x$fit)
	wh <- rep("f", length(newy))
	if(length(x$data) > 0) {
		dat <- x$data
		for(i in tv) {
			m <- x$xev[[i]]
			dat$x[, i] <- m[1 + round((dat$x[, i] - m[1])/(m[2] - m[
				1]))]
		}
		newx <- rbind(newx, dat$x)
		if(is.null(dat$y))
			newy <- c(newy, rep(NA, nrow(dat$x)))
		else {
			newy <- c(newy, dat$y/dat$sc)
			newy[is.na(newy)] <- 0
		}
		wh <- c(wh, rep("d", nrow(dat$x)))
	}
	if(length(tv) == 0) {
		newdat <- data.frame(newy, newx[, pv])
		names(newdat) <- c("y", paste("pv", 1:length(pv), sep = ""))
	}
	else {
		newdat <- data.frame(newx[, tv], newx[, pv], newy)
		names(newdat) <- c(paste("tv", 1:length(tv), sep = ""), paste(
			"pv", 1:length(pv), sep = ""), "y")
		for(i in 1:length(tv))
			newdat[, i] <- as.factor(signif(newdat[, i], 5))
	}
	loc.strip <- function(...)
	strip.default(..., strip.names = c(T, T), style = 1)
	if(length(pv) == 1) {
		clo <- cup <- numeric(0)
		if(length(x$se.fit) > 0) {
			if((!is.null(class(crit))) && (class(crit) == "kappa")
				)
				crit <- crit$crit.val
			cup <- x$trans((x$fit + crit * x$se.fit))
			clo <- x$trans((x$fit - crit * x$se.fit))
		}
		formula <- switch(1 + length(tv),
			y ~ pv1,
			y ~ pv1 | tv1,
			y ~ pv1 | tv1 * tv2,
			y ~ pv1 | tv1 * tv2 * tv3)
		pl <- xyplot(formula, xlab = pred.lab[pv], ylab = resp.lab, 
			main = main, type = "l", clo = clo, cup = cup, wh = wh, 
			panel = panel.xyplot.lf, data = newdat, strip = 
			loc.strip, ...)
	}
	if(length(pv) == 2) {
		formula <- switch(1 + length(tv),
			y ~ pv1 * pv2,
			y ~ pv1 * pv2 | tv1,
			y ~ pv1 * pv2 | tv1 * tv2,
			y ~ pv1 * pv2 | tv1 * tv2 * tv3)
		if(type == "contour")
			pl <- contourplot(formula, xlab = pred.lab[pv[1]], ylab
				 = pred.lab[pv[2]], main = main, data = newdat, 
				strip = loc.strip, ...)
		if(type == "level")
			pl <- levelplot(formula, xlab = pred.lab[pv[1]], ylab
				 = pred.lab[pv[2]], main = main, data = newdat, 
				strip = loc.strip, ...)
		if((type == "persp") | (type == "wireframe"))
			pl <- wireframe(formula, xlab = pred.lab[pv[1]], ylab
				 = pred.lab[pv[2]], zlab = resp.lab, data = 
				newdat, strip = loc.strip, ...)
	}
	if(length(tv) > 0)
		names(attr(pl$glist, "endpts")) <- attr(pl$glist, "names") <- 
			names(attr(pl$glist, "index")) <- pred.lab[tv]
	pl
}
"panel.xyplot.lf" <- 
function(x, y, subscripts, clo, cup, wh, type, ...)
{
	wh <- wh[subscripts]
	panel.xyplot(x[wh == "f"], y[wh == "f"], type = type, ...)
	if(length(clo) > 0) {
		panel.xyplot(x[wh == "f"], clo[subscripts][wh == "f"], type = 
			"l", lty = 2, ...)
		panel.xyplot(x[wh == "f"], cup[subscripts][wh == "f"], type = 
			"l", lty = 2, ...)
	}
	if(any(wh == "d")) {
		yy <- y[wh == "d"]
		if(any(is.na(yy)))
			rug(x[wh == "d"])
		else panel.xyplot(x[wh == "d"], yy)
	}
}
"plot.locfitpred" <- 
function(object, ...)
{
	if(object$dim == 1)
		plot.locfit.1d(object, ...)
	if(object$dim == 2)
		plot.locfit.2d(object, ...)
	if(object$dim >= 3)
		print(plot.locfit.3d(object, ...))
	invisible(NULL)
}
"panel.locfit" <- 
function(x, y, subscripts, z, xyz.labs, xyz.axes, xyz.mid, xyz.minmax, 
	xyz.range, col.regions, at, drape, contour, region, ...)
{
	if(!missing(z)) {
		z <- z[subscripts]
		fit <- locfit(z ~ x + y, ...)
		marg <- lfmarg(fit, m = 10)
		zp <- predict(fit, marg)
		if(!missing(contour))
			render.contour.trellis(marg[[1]], marg[[2]], zp, at = 
				at)
		else {
			loc.dat <- cbind(as.matrix(expand.grid(x = marg[[1]], y
				 = marg[[1]])), z = zp)
			print(xyz.minmax)
			print(xyz.range)
			print(xyz.mid)
			print(xyz.axes)
			render.3d.trellis(loc.dat, type = "wires", xyz.labs = 
				xyz.labs, xyz.axes = xyz.axes, xyz.mid = 
				xyz.mid, xyz.minmax = xyz.minmax, xyz.range = 
				xyz.range, col.regions = col.regions, at = at, 
				drape = drape)
		}
	}
	else {
		panel.xyplot(x, y)
		lines(locfit(y ~ x, ...))
	}
}
"lfmarg" <- 
function(xlim, m = 40)
{
	if(!is.numeric(xlim)) {
		d <- xlim$mi["d"]
		xlim <- xlim$box
	}
	else d <- length(m)
	marg <- vector("list", d)
	m <- rep(m, length.out = d)
	for(i in 1:d)
		marg[[i]] <- seq(xlim[i], xlim[i + d], length.out = m[i])
	marg
}
"lfeval" <- 
function(object)
object$eva
"lflim" <- 
function(limits, names, ret)
{
	d <- length(names)
	if(is.numeric(limits))
		ret <- limits
	else {
		z <- match(names, names(limits))
		for(i in 1:d)
			if(!is.na(z[i])) ret[c(i, i + d)] <- limits[[z[i]]]
	}
	as.numeric(ret)
}
"plot.eval" <- 
function(x, add = F, text = F, ...)
{
	d <- x$mi["d"]
	v <- matrix(x$eva$xev, nrow = d)
	ev <- x$mi["ev"]
	pv <- if(any(ev == c(1, 2))) as.logical(x$cell$s) else rep(F, ncol(v))
	if(!add) {
		plot(v[1,  ], v[2,  ], type = "n", xlab = x$vnames[1], ylab = x
			$vnames[2])
	}
	if(text)
		text(v[1,  ], v[2,  ], (1:x$nvc["nv"]) - 1)
	else {
		if(any(!pv))
			points(v[1, !pv], v[2, !pv], ...)
		if(any(pv))
			points(v[1, pv], v[2, pv], pch = "*", ...)
	}
	if(any(x$mi["ev"] == c(1, 2))) {
		zz <- .C("triterm",
			as.numeric(v),
			h = as.numeric(knots(x, what = "h", delete.pv = F)),
			as.integer(x$cell$ce),
			lo = as.integer(x$cell$lo),
			hi = as.integer(x$cell$hi),
			as.numeric(x$eva$scale),
			as.integer(x$nvc),
			as.integer(x$mi),
			as.numeric(x$dp),
			nt = integer(1),
			term = integer(600))
		ce <- zz$term + 1
	}
	else ce <- x$cell$ce + 1
	if(any(x$mi["ev"] == c(1, 5, 7))) {
		vc <- 2^d
		ce <- matrix(ce, nrow = vc)
		segments(v[1, ce[1,  ]], v[2, ce[1,  ]], v[1, ce[2,  ]], v[2, 
			ce[2,  ]], ...)
		segments(v[1, ce[1,  ]], v[2, ce[1,  ]], v[1, ce[3,  ]], v[2, 
			ce[3,  ]], ...)
		segments(v[1, ce[2,  ]], v[2, ce[2,  ]], v[1, ce[4,  ]], v[2, 
			ce[4,  ]], ...)
		segments(v[1, ce[3,  ]], v[2, ce[3,  ]], v[1, ce[4,  ]], v[2, 
			ce[4,  ]], ...)
	}
	if(any(x$mi["ev"] == c(2, 8))) {
		vc <- d + 1
		m <- matrix(ce, nrow = 3)
		segments(v[1, m[1,  ]], v[2, m[1,  ]], v[1, m[2,  ]], v[2, m[2, 
			]], ...)
		segments(v[1, m[1,  ]], v[2, m[1,  ]], v[1, m[3,  ]], v[2, m[3, 
			]], ...)
		segments(v[1, m[2,  ]], v[2, m[2,  ]], v[1, m[3,  ]], v[2, m[3, 
			]], ...)
	}
	invisible(NULL)
}
"rv" <- 
function(fit)
fit$dp["rv"]
"rv<-" <- 
function(fit, value)
{
	fit$dp["rv"] <- value
	fit
}
"summary.locfit" <- 
function(object, ...)
{
	mi <- object$mi
	fam <- c("Density Estimation", "Poisson process rate estimation", 
		"Hazard Rate Estimation", "Local Regression", 
		"Local Likelihood - Binomial", "Local Likelihood - Poisson", 
		"Local Likelihood - Gamma", "Local Likelihood - Geometric", 
		"Local Robust Regression")[mi["tg"]]
	estr <- c("Rectangular Tree", "Triangulation", "Data", 
		"Rectangular Grid", "k-d tree", "k-d centres", 
		"Cross Validation", "User-provided")[mi["ev"]]
	ret <- list(call = object$call, etype = fam, n = mi["n"], d = mi["d"], 
		estr = estr, nv = object$nvc["nv"], deg = mi["deg"], dp = 
		object$dp, vnames = object$vnames)
	class(ret) <- "summary.locfit"
	ret
}
"summary.locfitpred" <- 
function(object, ...)
object$trans(object$fit)
"summary.locfit.1d" <- 
function(object, ...)
object$trans(object$fit)
"summary.locfit.2d" <- 
function(object, ...)
object$trans(object$fit)
"summary.locfit.2d" <- 
function(object, ...)
object$trans(object$fit)
"print.summary.locfit" <- 
function(x)
{
	cat("Estimation type:", x$etype, "\n")
	cat("\nCall:\n")
	print(x$call)
	cat("\nNumber of data points: ", x$n, "\n")
	cat("Independent variables: ", x$vnames, "\n")
	cat("Evaluation structure:", x$estr, "\n")
	cat("Number of evaluation points: ", x$nv, "\n")
	cat("Degree of fit: ", x$deg, "\n")
	cat("Fitted Degrees of Freedom: ", round(x$dp["t1"], 3), "\n")
	invisible(x)
}
"hatmatrix" <- 
function(formula, dc = T, ...)
{
	m <- match.call()
	m$geth <- 1
	m[[1]] <- as.name("locfit")
	eval(m, sys.frame(sys.parent()))
}
"regband" <- 
function(formula, what = c("CP", "GCV", "GKK", "RSW"), deg = 1, ...)
{
	m <- match.call()
	m$geth <- 3
	m$deg <- c(deg, 4)
	m$what <- NULL
	m$deriv <- match(what, c("CP", "GCV", "GKK", "RSW"))
	m[[1]] <- as.name("locfit")
	z <- eval(m, sys.frame(sys.parent()))
	names(z) <- what
	z
}
"kdeb" <- 
function(x, h0 = 0.01 * sd, h1 = sd, what = c("AIC", "LCV", "LSCV", "BCV", 
	"SJPI", "GKK"), kern = "gauss", gf = 2.5)
{
	n <- length(x)
	sd <- sqrt(var(x))
	z <- .C("kdeb",
		x = as.numeric(x),
		mi = as.integer(n),
		band = numeric(length(what)),
		ind = integer(n),
		h0 = as.numeric(gf * h0),
		h1 = as.numeric(gf * h1),
		meth = as.integer(match(what, c("AIC", "LCV", "LSCV", "BCV", 
			"SJPI", "GKK"))),
		nmeth = as.integer(length(what)),
		kern = pmatch(kern, c("rect", "epan", "bisq", "tcub", "trwt", 
			"gauss")))
	band <- z$band
	names(band) <- what
	band
}
"knots" <- 
function(x, tr, what = c("x", "coef", "h", "nlx"), delete.pv = T, ...)
{
	nv <- x$nvc["nv"]
	d <- x$mi["d"]
	p <- x$mi["p"]
	z <- 0:(nv - 1)
	ret <- matrix(0, nrow = nv, ncol = 1)
	rname <- character(0)
	if(missing(tr))
		tr <- x$trans
	coef <- x$eva$coef
	for(wh in what) {
		if(wh == "x") {
			ret <- cbind(ret, matrix(x$eva$xev, ncol = d, byrow = T
				))
			rname <- c(rname, x$vnames)
		}
		if(wh == "coef") {
			d0 <- coef[, 1]
			d0[d0 == 0.1278433] <- NA
			ret <- cbind(ret, tr(d0))
			rname <- c(rname, "mu hat")
		}
		if(wh == "f1") {
			ret <- cbind(ret, coef[, 1 + (1:d)])
			rname <- c(rname, paste("d.", 1:d, sep = ""))
		}
		if(wh == "f2") {
			dd <- (d * (d + 1))/2
			d2 <- coef[, 1 + d + 1:dd]
			u <- matrix(0, nrow = d, ncol = d)
			cu <- col(u)
			ru <- row(u)
			ret <- cbind(ret, matrix(d2, ncol = dd))
			rname <- c(rname, paste("d.", cu[cu <= ru], ".", ru[cu <= 
				ru], sep = ""))
		}
		if(wh == "coefa") {
			ret <- cbind(ret, coef[, 1:p])
			rname <- c(rname, paste("c.", 1:p, sep = ""))
		}
		if(wh == "nlx") {
			ret <- cbind(ret, coef[, p + 1])
			rname <- c(rname, "||l(x)||")
		}
		if(wh == "nlxa") {
			ret <- cbind(ret, coef[, p + 1 + (0:d)])
			rname <- c(rname, "nlx", paste("nlx-d", 1:d, sep = ""))
		}
		if(wh == "se") {
			ret <- cbind(ret, sqrt(x$dp["rv"]) * coef[, p + 1])
			rname <- c(rname, "StdErr")
		}
		if(wh == "infl") {
			z <- coef[, 2 * p + d + 1]
			ret <- cbind(ret, z * z)
			rname <- c(rname, "Influence")
		}
		if(wh == "infla") {
			ret <- cbind(ret, coef[, 2 * p + d + 1 + (0:d)])
			rname <- c(rname, "inf^.5", paste("inf-d", 1:d, sep = 
				""))
		}
		if(wh == "lik") {
			ret <- cbind(ret, coef[, 2 * p + 2 * d + 2:4])
			rname <- c(rname, c("LocLike", "fit.df", "res.df"))
		}
		if(wh == "h") {
			ret <- cbind(ret, coef[, 2 * p + 2 * d + 5])
			rname <- c(rname, "h")
		}
		if(wh == "deg") {
			ret <- cbind(ret, coef[, 2 * p + 2 * d + 6])
			rname <- c(rname, "deg")
		}
	}
	ret <- as.matrix(ret[, -1])
	if(nv == 1)
		ret <- t(ret)
	dimnames(ret) <- list(NULL, rname)
	if((delete.pv) && (any(x$mi["ev"] == c(1, 2))))
		ret <- ret[!as.logical(x$cell$s),  ]
	ret
}
"locfit.matrix" <- 
function(fit, data)
{
	m <- fit$call
	if(m[[1]] == "locfit.raw") {
		x <- as.matrix(eval(m$x, fit$frame))
		y <- if(is.null(m$y)) 0 else eval(m$y, fit$frame)
		w <- if(is.null(m$weights)) 1 else eval(m$weights, fit$frame)
		ce <- eval(m$cens, fit$frame)
		if(is.null(ce))
			ce <- 0
		base <- if(is.null(m$base)) 0 else eval(m$base, fit$frame)
	}
	else {
		Terms <- terms(as.formula(m$formula))
		attr(Terms, "intercept") <- 0
		m[[1]] <- as.name("model.frame")
		z <- pmatch(names(m), c("formula", "data", "weights", "cens", 
			"base", "subset"))
		for(i in length(z):2)
			if(is.na(z[i])) m[[i]] <- NULL
		frm <- eval(m, fit$frame)
		vnames <- as.character(attributes(Terms)$variables)[-1]
		if(attr(Terms, "response")) {
			y <- model.extract(frm, response)
			vnames <- vnames[-1]
		}
		else y <- 0
		x <- as.matrix(frm[, vnames])
		if(any(names(m) == "weights"))
			w <- model.extract(frm, weights)
		else w <- 1
		if(any(names(m) == "cens"))
			ce <- model.extract(frm, cens)
		else ce <- 0
		if(any(names(m) == "base"))
			base <- model.extract(frm, base)
		else base <- 0
	}
	sc <- if(any(fit$mi["tg"] == 5:8)) w else 1
	list(x = x, y = y, w = w, sc = sc, ce = ce, base = base)
}
"kappa0" <- 
function(formula, dc = T, cov = 0.95, ...)
{
	if(class(formula) == "locfit") {
		m <- formula$call
	}
	else {
		m <- match.call()
		m$cov <- NULL
	}
	m$dc <- T
	m$geth <- 2
	m[[1]] <- as.name("locfit")
	z <- eval(m, sys.frame(sys.parent()))
	z$cov <- cov
	z$crit.val <- .C("compcv",
		k0 = z$const,
		d = as.integer(z$d),
		cov = as.numeric(cov),
		x = numeric(1))$x
	class(z) <- "kappa"
	z
}
"gcv" <- 
function(...)
{
	m <- match.call()
	m[[1]] <- as.name("locfit")
	fit <- eval(m, sys.frame(sys.parent()))
	z <- fit$dp[c("lk", "t0", "t1")]
	n <- fit$mi["n"]
	z <- c(z, (-2 * n * z[1])/(n - z[2])^2)
	names(z) <- c("lik", "infl", "vari", "gcv")
	z
}
"aic" <- 
function(..., pen = 2)
{
	m <- match.call()
	m[[1]] <- as.name("locfit")
	m$pen <- NULL
	dp <- eval(m, sys.frame(sys.parent()))$dp
	z <- dp[c("lk", "t0", "t1")]
	z <- c(z, -2 * z[1] + pen * z[2])
	names(z) <- c("lik", "infl", "vari", "aic")
	z
}
"dv" <- 
function(z)
{
	n <- length(z)
	e <- z[ - c(1, n)] - (z[ - c(1, 2)] + z[ - c(n - 1, n)])/2
	return((2 * sum(e * e))/(3 * (n - 2)))
}
"expit" <- 
function(x)
{
	y <- numeric(length(x))
	y[x > 0] <- 1/(1 + exp( - x[x > 0]))
	z <- exp(x[x <= 0])
	y[x <= 0] <- z/(1 + z)
	y
}
"fw" <- 
function(x, k = 1)
{
	f <- rep(0, length(x))
	if(k == 1) {
		tj <- c(0.10000000000000001, 0.13, 0.14999999999999999, 
			0.23000000000000001, 0.25, 0.40000000000000002, 0.44, 
			0.65000000000000002, 0.76000000000000001, 
			0.78000000000000003, 0.81000000000000005)
		hj <- c(4, -5, 3, -4, 5, -4.2000000000000002, 
			2.1000000000000001, 4.2999999999999998, 
			-3.1000000000000001, 2.1000000000000001, 
			-4.2000000000000002)
		for(j in 1:length(tj))
			f <- f + 3.6000000000000001 * (hj[j] * (x > tj[j]))
	}
	if(k == 2) {
		tj <- c(0.10000000000000001, 0.13, 0.14999999999999999, 
			0.23000000000000001, 0.25, 0.40000000000000002, 0.44, 
			0.65000000000000002, 0.76000000000000001, 
			0.78000000000000003, 0.81000000000000005)
		hj <- c(4, 5, 3, 4, 5, 4.2000000000000002, 2.1000000000000001, 
			2.1000000000000001, 4.2999999999999998, 
			-3.1000000000000001, 2.1000000000000001, 
			-4.2000000000000002)
		for(j in 1:length(tj))
			f <- f + 3.6000000000000001 * (hj[j] * (x > tj[j]))
	}
	if(k == 2) {
		tj <- c(0.10000000000000001, 0.13, 0.14999999999999999, 
			0.23000000000000001, 0.25, 0.40000000000000002, 0.44, 
			0.65000000000000002, 0.76000000000000001, 
			0.78000000000000003, 0.81000000000000005)
		hj <- c(4, 5, 3, 4, 5, 4.2000000000000002, 2.1000000000000001, 
			4.2999999999999998, 3.1000000000000001, 
			5.0999999999999996, 4.2000000000000002)
		wj <- c(5, 5, 6, 10, 10, 30, 10, 10, 5, 8, 5)/1000
		for(j in 1:length(tj))
			f <- f + (10 * hj[j])/(1 + abs((x - tj[j])/wj[j]))^4
	}
	if(k == 3) {
		f <- 10 * sin(4 * pi * x) - 5 * ((x > 0.29999999999999999) & (x <
			0.71999999999999997))
	}
	if(k == 4) {
		f <- 20 * sqrt(x * (1 - x)) * sin((2 * pi * 1.05)/(x + 
			0.050000000000000003))
	}
	f
}
"locfit.robust" <- 
function(x, y, weights, ..., iter = 3)
{
	n <- length(y)
	lfr.wt <- rep(1, n)
	for(i in 0:iter) {
		fit <- locfit.raw(x, y, weights = lfr.wt, ...)
		res <- residuals(fit)
		s <- median(abs(res))
		lfr.wt <- pmax(1 - (res/(6 * s))^2, 0)^2
	}
	fit$frame <- sys.frame(sys.parent())
	lfr.wt <<- lfr.wt
	fit
}
"locfit.censor" <- 
function(x, y, cens, ..., iter = 3, km = F)
{
	lfc.y <- y
	cens <- as.logical(cens)
	for(i in 0:iter) {
		fit <- locfit.raw(x, lfc.y, ...)
		fh <- fitted(fit)
		if(km) {
			sr <- y - fh
			lfc.y <- fh + km.mrl(sr, cens)
		}
		else {
			rdf <- sum(1 - cens) - 2 * fit$dp["t0"] + fit$dp["t1"]
			sigma <- sqrt(sum((y - fh) * (lfc.y - fh))/rdf)
			sr <- (y - fh)/sigma
			lfc.y <- fh + (sigma * dnorm(sr))/pnorm( - sr)
		}
		lfc.y[!cens] <- y[!cens]
	}
	fit$frame <- sys.frame(sys.parent())
	lfc.y <<- lfc.y
	fit
}
"locfit.quasi" <- 
function(x, y, weights, ..., iter = 3, var = function(mean)
abs(mean))
{
	n <- length(y)
	w0 <- lfq.wt <- if(missing(weights)) rep(1, n) else weights
	for(i in 0:iter) {
		fit <- locfit.raw(x, y, weights = lfq.wt, ...)
		fh <- fitted(fit)
		lfq.wt <- w0/var(fh)
	}
	fit$frame <- sys.frame(sys.parent())
	lfq.wt <<- lfq.wt
	fit
}
"locfit.add" <- 
function(x1, x2, y, iter = 5)
{
	base <- 0
	for(i in 1:iter) {
		assign("loc.dat", data.frame(x1 = x1, y = y, base = base), 
			frame = 1)
		fit1 <- locfit(y ~ x1, base = base, data = loc.dat)
		base <- fitted(fit1, data = loc.dat)
		assign("loc.dat", data.frame(x2 = x2, y = y, base = base), 
			frame = 1)
		fit2 <- locfit(y ~ x2, base = base, data = loc.dat)
		base <- fitted(fit2, data = loc.dat)
		base <- base - base[1]
		print(fit1$dp["lk"])
	}
	list(f1 = fit1, f2 = fit2)
}
"plotbyfactor" <- 
function(x, y, f, data, col = 1:10, pch = "O", add = F, lg, xlab = deparse(
	substitute(x)), ylab = deparse(substitute(y)), ...)
{
	if(!missing(data)) {
		x <- eval(substitute(x), data)
		y <- eval(substitute(y), data)
		f <- eval(substitute(f), data)
	}
	f <- as.factor(f)
	if(!add)
		plot(x, y, type = "n", xlab = xlab, ylab = ylab, ...)
	lv <- levels(f)
	col <- rep(col, length.out = length(lv))
	pch <- rep(pch, length.out = length(lv))
	for(i in 1:length(lv)) {
		ss <- f == lv[i]
		if(any(ss))
			points(x[ss], y[ss], col = col[i], pch = pch[i])
	}
	if(!missing(lg))
		legend(lg[1], lg[2], legend = levels(f), col = col, pch = paste(
			pch, collapse = ""))
}