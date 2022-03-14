# Internal function that gets metadata from outputs
get_metadata <- function(outputs) {

    # If metadata is stored with outputs, return it
    if (!is.null(outputs$meta)) return(outputs$meta)

    # Get all names recursively
    all_leaf_names <- function(x) {
        if (!is.list(x)) return(character(0))
        ll <- sapply(x, is.list)
        c(names(x)[!ll], as.character(unlist(lapply(x, all_leaf_names))))
    }

    nm <- unique(all_leaf_names(outputs))

    list(parameters=lapply(nm, function(x) list(name=x)))
}

#' Create a `data.frame` of from outputs and metadata
#'
#' This can be viewed as the first step in creating a nice-looking HTML table
#' of model parameters. It combines the "raw" model outputs with metadata and
#' produces and `data.frame`, conceived as an intermediate between the raw
#' outputs and formatted table, but may also be useful in its own right. The
#' decoupling of raw outputs from the final table is viewed as essential for
#' flexibility.
#'
#' @param outputs A `list` of outputs from fitting the model (see Details).
#' @param meta A `list` of metadata (see Details).
#'
#' @details One of the key features of the approach taken in this package is
#' that it decouples the "raw" outputs of the model from the presentation of
#' results. A metadata description of the desired presentation of results is
#' what links the two.  This allows, for example, parameters to be presented in
#' a different order, or on a different scale, than they were specified in the
#' model. Hence, it provides more flexibility and control over the presentation
#' than other approaches.
#'
#' `outputs` is a named `list`, with the following elements:
#' - `est`:       estimated values (i.e., point estimates)
#' - `se`:        standard errors
#' - `fixed`:     designates parameters that were fixed rather than estimated
#' - `shrinkage`: for random effects, the estimated percent shrinkage
#'
#' `est`, `se` and `fixed` have essentially the same structure. They can be
#' either flat named lists, or more structured named lists containing the following elements:
#'
#' - `th`     : named list (or vector) of fixed effects
#' - `om`     : named list (or vector) of individual-level random effects
#' expressed as standard deviations
#' - `om_cov` : individual-level random effects expressed as a
#' variance-covariance matrix
#' - `om_cor` : individual-level random effects expressed as a matrix of
#' correlations (off-diagonal elements) and standard deviations (diagonal
#' elements)
#' - `sg`     : named list (or vector) of observation-level random effects
#' expressed as standard deviations
#'
#' `meta` is a `list`. Each element of `meta` is a named (sub)list
#' representing a parameter.  Each parameter is described by a series of
#' attributes (not R `attributes`, but named list items). Of these, the only
#' one that is required is `name`, which must match the name of the parameter
#' used in `outputs` as it is used to make that association. The optional attributes include:
#'
#' - `label`:     A descriptive label.
#' - `units`:     Units, if applicable.
#' - `type`:      Parameters can be grouped into sections by type. The standard types are:
#'   - `Structural`:      Structural model parameters
#'   - `CovariateEffect`: Parameters that relate covariates to structural parameters
#'   - `IIV`:             Inter-individual (i.e., between-subject) variability
#'   - `IOV`:             Inter-occasion variability
#'   - `RUV`:             Residual unexplained variability
#' - `trans`:     Parameters can be presented on a (back)transformed scale (e.g.,
#'                antilog). Importantly, transformation are also applied to
#'                standard errors (by "propagation of errors", also known as
#'                the delta method) to preserve (asymptotic) correctness, and
#'                to the endpoints of confidence intervals (note: this
#'                typically leads to non-symmetric intervals). Only a small set
#'                of transformations are currently recognized and supported,
#'                which include:
#'   - `identity`:  no transformation
#'   - `%`,         percent-scale
#'   - `exp`:       antilog
#'   - `ilogit`:    inverse-logit
#'   - `CV%`:       intended specifically for IIV parameters, where the
#'                  associated structural parameter is log-normally
#'                  distributed, transforms the standard deviation \eqn{\omega}
#'                  to percent coefficient of variation by the formula
#'                  \eqn{100\times\sqrt{\exp(\omega^2)-1}}
#'   - `SD (CV%)`:  similar to the above, but the parameter remains on its original
#'                  scale (i.e., standard deviation) with the percent coefficient
#'                  of variation displayed in parentheses next to it (does not
#'                  affect standard errors or confidence intervals).
#'
#' @return A `data.frame` with a row for each parameter, and the following columns:
#'
#' - `name`:      name of the parameter (`character`)
#' - `fixed`:     fixed or estimated? (`logical`)
#' - `est`:       estimated value (`numeric`)
#' - `se`:        standard error (`numeric`)
#' - `rse`:       percent relative standard error (`numeric`)
#' - `lci95`:     lower bound of 95% confidence interval (`numeric`)
#' - `uci95`:     upper bound of 95% confidence interval (`numeric`)
#' - `pval`:      p-value for test of null hypothesis that value is zero (`numeric`)
#' - `shrinkage`: percent shrinkage if applicable (`numeric`)
#'
#' Other attributes from `meta` will also be preserved as columns. The order of
#' the rows is determined by the order of the parameters in `meta` (the order
#' in `outputs` is irrelevant).
#'
#' @seealso [pmxpartab]
#'
#' @examples
#' outputs <- list(
#'   est = list(
#'     th = list(CL = 0.482334, VC = 0.0592686),
#'     om = list(nCL = 0.315414, nVC = 0.536025),
#'     sg = list(ERRP = 0.0508497)),
#'   se = list(
#'     th = list(CL = 0.0138646, VC = 0.00555121),
#'     om = list(nCL = 0.0188891, nVC = 0.0900352),
#'     sg = list(ERRP = 0.00182851)),
#'   fixed = list(
#'     th = list(CL = FALSE, VC = FALSE),
#'     om = list(nCL = FALSE, nVC = FALSE),
#'     sg = list(ERRP = FALSE)),
#'   shrinkage = list(nCL = 9.54556, nVC = 47.8771))
#' 
#' meta <- list(
#'   parameters = list(
#'     list(name="CL", label="Clearance", units="L/h", type="Structural"),
#'     list(name="VC", label="Volume", units="L", type="Structural", trans="exp"),
#'     list(name="nCL", label="On Clearance", type="IIV", trans="SD (CV%)"),
#'     list(name="nVC", label="On Volume", type="IIV"),
#'     list(name="ERRP", label="Proportional Error", units="%", type="RUV", trans="%")))
#' 
#' pmxparframe(outputs, meta)
#' @importFrom stats pnorm
#' @export
pmxparframe <- function(outputs, meta=get_metadata(outputs)) {
    param <- meta$parameters
    z <- data.table::rbindlist(param, fill=T)

    # Legacy
    if (!is.null(outputs$th))      outputs$est$th <- outputs$th
    if (!is.null(outputs$sg))      outputs$est$sg <- outputs$sg
    if (!is.null(outputs$om))      outputs$est$om <- outputs$om

    if (is.null(outputs$est))      outputs$est      <- list()
    if (is.null(outputs$est$th))   outputs$est$th   <- list()
    if (is.null(outputs$est$om))   outputs$est$om   <- list()
    if (is.null(outputs$est$sg))   outputs$est$sg   <- list()
    if (is.null(outputs$se))       outputs$se       <- list()
    if (is.null(outputs$se$th))    outputs$se$th    <- list()
    if (is.null(outputs$se$om))    outputs$se$om    <- list()
    if (is.null(outputs$se$sg))    outputs$se$sg    <- list()
    if (is.null(outputs$fixed))    outputs$fixed    <- list()
    if (is.null(outputs$fixed$th)) outputs$fixed$th <- list()
    if (is.null(outputs$fixed$om)) outputs$fixed$om <- list()
    if (is.null(outputs$fixed$sg)) outputs$fixed$sg <- list()

    nm <- setdiff(names(outputs$est), c("th", "om", "om_cov", "om_cor", "sg", "all"))
    outputs$est$all <- c(outputs$est$all, outputs$est[nm])

    nm <- setdiff(names(outputs$se), c("th", "om", "om_cov", "om_cor", "sg", "all"))
    outputs$se$all <- c(outputs$se$all, outputs$se[nm])

    nm <- setdiff(names(outputs$fixed), c("th", "om", "om_cov", "om_cor", "sg", "all"))
    outputs$fixed$all <- c(outputs$fixed$all, outputs$fixed[nm])

    outputs$est$all   <- with(outputs$est, c(all, th, om, sg))
    outputs$se$all    <- with(outputs$se, c(all, th, om, sg))
    outputs$fixed$all <- with(outputs$fixed, c(all, th, om, sg))

    z$fixed     <- as.logical(NA)
    z$est       <- as.numeric(NA)
    z$se        <- as.numeric(NA)
    z$rse       <- as.numeric(NA)
    z$lci95     <- as.numeric(NA)
    z$uci95     <- as.numeric(NA)
    z$pval      <- as.numeric(NA)
    if (!is.null(outputs$shrinkage)) {
        z$shrinkage <- as.numeric(NA)
    }

    have.bootstrap <- !is.null(outputs$bootstrap)
    if (have.bootstrap) {
        z$boot.median <- as.numeric(NA)
        z$boot.lci    <- as.numeric(NA)
        z$boot.uci    <- as.numeric(NA)
    }

    for (i in 1:nrow(z)) {

        `%||%` <- function(x, y) { if (is.null(x) || is.na(x)) y else x }
        name <- z$name[i] %||% NA
        trans <- z$trans[i] %||% NA

        if (have.bootstrap) {
            boot.median <- NA
            boot.lci <- NA
            boot.uci <- NA
        }

        # Check parameter type
        if (name %in% names(outputs$est$all)) {
            est <- outputs$est$all[[name]]
            se <- outputs$se$all[[name]]
            fixed <- outputs$fixed$all[[name]]
            if (have.bootstrap && name %in% names(outputs$bootstrap$median)) {
                boot.median <- outputs$bootstrap$median[[name]]
                boot.lci <- outputs$bootstrap$ci95[[name]][1]
                boot.uci <- outputs$bootstrap$ci95[[name]][2]
            }
        } else if ((re <- regexec("^om(?:_cor)?\\((\\w+),(\\w+)\\)$", name, perl=T))[[1]][1] > 0) {
            a1 <- regmatches(name, re)[[1]][2]
            a2 <- regmatches(name, re)[[1]][3]
            if (is.null(dimnames(outputs$est$om_cor)) && !is.null(names(outputs$se$om))) {
                dimnames(outputs$est$om_cor) <- list(names(outputs$est$om), names(outputs$est$om))
            }
            if (is.null(dimnames(outputs$se$om_cor)) && !is.null(names(outputs$se$om))) {
                dimnames(outputs$se$om_cor) <- list(names(outputs$se$om), names(outputs$se$om))
            }
            if (a1 %in% dimnames(outputs$est$om_cor)[[1]] && a2 %in% dimnames(outputs$est$om_cor)[[2]]) {
                est <- outputs$est$om_cor[a1,a2]
                se <- outputs$se$om_cor[a1,a2]
                fixed <- outputs$fixed$om_cor[a1,a2]
                if (isTRUE(fixed)) {
                    se <- NA
                }
                if (have.bootstrap) {
                    boot.re <- paste0("OMEGA\\D", max(a1, a2), "\\D", min(a1, a2), "($|\\D)")
                    boot.name <- grep(boot.re, names(outputs$bootstrap$median), ignore.case=T, value=T)
                    if (length(boot.name) == 1) {
                        boot.median <- outputs$bootstrap$median[[boot.name]]
                        boot.lci <- outputs$bootstrap$ci95[[boot.name]][1]
                        boot.uci <- outputs$bootstrap$ci95[[boot.name]][2]
                    }
                }
            } else {
                # Not found
                next
            }
        } else if ((re <- regexec("^om_cov\\((\\w+),(\\w+)\\)$", name, perl=T))[[1]][1] > 0) {
            a1 <- regmatches(name, re)[[1]][2]
            a2 <- regmatches(name, re)[[1]][3]
            if (is.null(dimnames(outputs$est$om_cov))) {
                dimnames(outputs$est$om_cov) <- list(names(outputs$est$om), names(outputs$est$om))
            }
            if (is.null(dimnames(outputs$se$om_cov)) && !is.null(names(outputs$se$om))) {
                dimnames(outputs$se$om_cov) <- list(names(outputs$se$om), names(outputs$se$om))
            }
            if (a1 %in% dimnames(outputs$est$om_cov)[[1]] && a2 %in% dimnames(outputs$est$om_cov)[[2]]) {
                est <- outputs$est$om_cov[a1,a2]
                se <- outputs$se$om_cov[a1,a2]
                fixed <- outputs$fixed$om_cov[a1,a2]
                if (isTRUE(fixed)) {
                    se <- NA
                }
                if (have.bootstrap) {
                    boot.re <- paste0("OMEGA\\D", max(a1, a2), "\\D", min(a1, a2), "($|\\D)")
                    boot.name <- grep(boot.re, names(outputs$bootstrap$median), ignore.case=T, value=T)
                    if (length(boot.name) == 1) {
                        boot.median <- outputs$bootstrap$median[[boot.name]]
                        boot.lci <- outputs$bootstrap$ci95[[boot.name]][1]
                        boot.uci <- outputs$bootstrap$ci95[[boot.name]][2]
                    }
                }
            } else {
                # Not found
                next
            }
        } else {
            # Not found
            next
        }

        if (is.null(fixed)) {
            fixed <- FALSE
        }

        if (isTRUE(fixed) || is.null(se) || is.na(se)) {
            se <- NA
            rse <- NA
            ci95 <- c(NA, NA)
            pval <- NA
        } else {
            rse <- 100*se/abs(est)
            ci95 <- est + c(-1,1)*1.96*se
            pval <- 2*(1 - pnorm(abs(est/se)))
        }

        # Check transformation
        if (!is.na(trans) && trans == "%") {
            est <- 100*est
            if (!fixed) {
                se  <- 100*se
                ci95  <- 100*ci95
            }
            if (have.bootstrap) {
                boot.median <- 100*boot.median
                boot.lci <- 100*boot.lci
                boot.uci <- 100*boot.uci
            }
        } else if (!is.na(trans) && trans == "exp") {
            est <- exp(est)
            if (!fixed) {
                rse <- 100*se
                se  <- (rse/100)*est
                ci95  <- exp(ci95)
            }
            if (have.bootstrap) {
                boot.median <- exp(boot.median)
                boot.lci <- exp(boot.lci)
                boot.uci <- exp(boot.uci)
            }
        } else if (!is.na(trans) && trans == "ilogit") {
            ilogit <- function(x) { 1 / (1 + exp(-x)) }
            est <- ilogit(est)
            if (!fixed) {
                rse <- 100*se*(1 - est)
                se  <- (rse/100)*est
                ci95  <- ilogit(ci95)
            }
            if (have.bootstrap) {
                boot.median <- ilogit(boot.median)
                boot.lci <- ilogit(boot.lci)
                boot.uci <- ilogit(boot.uci)
            }
        } else if (!is.na(trans) && trans == "CV%") {
            g <- function(x) { 100*sqrt(exp(x^2) - 1) }
            dg <- function(x) { 100*0.5*(1/sqrt(exp(x^2) - 1))*exp(x^2)*2*x }
            x <- est
            est <- g(x)
            if (!fixed) {
                se  <- se*dg(x)
                rse <- 100*se/abs(est)
                ci95  <- g(ci95)
            }
            if (have.bootstrap) {
                boot.median <- g(boot.median)
                boot.lci <- g(boot.lci)
                boot.uci <- g(boot.uci)
            }
        }

        z$fixed[i] <- fixed
        z$est[i]   <- est
        z$se[i]    <- se
        z$rse[i]   <- rse
        z$lci95[i] <- ci95[1]
        z$uci95[i] <- ci95[2]
        z$pval[i]  <- pval
        if (have.bootstrap) {
            z$boot.median[i] <- boot.median
            z$boot.lci[i] <- boot.lci
            z$boot.uci[i] <- boot.uci
        }

        if (name %in% names(outputs$shrinkage)) {
            z$shrinkage[i] <- outputs$shrinkage[[name]]
        }
    }
    z <- z[!is.na(z$est),, drop=FALSE]
    as.data.frame(z)
}

# Internal function to help format numbers
p <- function(x, digits=3, flag="", round.integers=FALSE){
    if (!is.numeric(x)) {
        return(x)
    }
    prefix <- ifelse(flag=="+" & x > 0, "+", "")
    x <- ifelse(!is.na(x) & x == Inf, "\U{221e}",
         ifelse(!is.na(x) & x == -Inf, "-\U{221e}",
        table1::signif_pad(x, digits=digits, round.integers=round.integers)))
    paste0(prefix, x)
}

# Internal function that produces a table section heading
partab_section <- function(label, ncolumns) {
    paste0(c('<tr>',
        paste0(sprintf('<td class="partabsectionheading">%s</td>', c(label, rep("", ncolumns-1))), collapse='\n'),
        '</tr>'), collapse='\n')
}

# Internal function that produces a single table row
partab_row <- function(
    name,
    label          = NULL,
    units          = NULL,
    type           = NULL,
    trans          = NULL,
    expression     = NULL,
    relatedTo      = NULL,
    superscript    = NULL,
    fixed          = NULL,
    est            = NULL,
    se             = NULL,
    rse            = NULL,
    pval           = NULL,
    lci95          = NULL,
    uci95          = NULL,
    boot.median    = NULL,
    boot.lci       = NULL,
    boot.uci       = NULL,
    shrinkage      = NULL,
    merge.units    = TRUE,
    na             = "-",
    digits         = 3,
    indent         = TRUE,
    have.bootstrap = !is.null(boot.median),
    columns        = c(est="Estimate", rse="RSE%", ci95="95% CI", shrinkage="Shrinkage"),
    ...) {

    # Check for superscript
    if (is.null(superscript) || is.na(superscript)) {
        superscript <- ""
    } else {
        superscript <- paste0("<sup>", superscript, "</sup>")
    }

    # Check for label
    if (is.null(label) || is.na(label)) {
        label <- name
    }

    # Check for units
    if (!is.null(units) && !is.na(units)) {
        if (merge.units) {
            label <- sprintf("%s (%s)", label, units)
        }
    }

    if (!is.null(trans) && !is.na(trans) && trans == "SD (CV%)") {
        g <- function(x) { 100*sqrt(exp(x^2) - 1) }
        x <- est
        est <- sprintf("%s (%s%%)", p(x, digits), p(g(x), digits))
    } else {
        est <- p(est, digits)
    }

    est <- paste0(est, superscript)

    if (isTRUE(fixed)) {
        est <- sprintf('%s Fixed', est)
    }

    if (is.null(se) || is.na(se)) {
        se <- na
    } else {
        se <- p(se, digits)
    }

    if (is.null(rse) || is.na(rse)) {
        rse <- na
    } else {
        rse <- p(rse, digits)
    }

    if (is.null(lci95) || is.null(uci95)) {
        ci95 <- na
    } else if (is.na(lci95) || is.na(uci95)) {
        ci95 <- na
    } else {
        ci95 <- sprintf('%s &ndash; %s', p(lci95, digits), p(uci95, digits))
    }

    if (is.null(pval) || is.na(pval)) {
        pval <- na
    } else {
        pval <- fpval(pval, digits=digits, html=TRUE)
    }

    if (have.bootstrap) {
        if (is.na(boot.median)) {
            boot.median <- na
            boot.ci95 <- na
        } else {
            boot.median <- p(boot.median, digits)
            boot.ci95 <- sprintf('%s &ndash; %s', p(boot.lci, digits), p(boot.uci, digits))
        }
    } else {
        boot.ci95 <- NULL
    }

    if (!is.null(shrinkage)) {
        if (is.na(shrinkage)) {
            shrinkage <- na
        } else {
            shrinkage <- sprintf("%s%%", p(shrinkage, digits))
        }
    }
    all <- c(est=est, se=se, rse=rse, pval=pval, ci95=ci95, boot.median=boot.median, boot.ci95=boot.ci95, shrinkage=shrinkage)
    args <- list(...)
    if (!is.null(names(args))) {
        all <- c(all, args[!(names(args) %in% c("", names(all)))])
    }
    paste0(c('<tr>',
        sprintf('<td class="%s">%s</td>', ifelse(isTRUE(indent), "partablabelindent", "partablabelnoindent"), label),
        paste0(sprintf('<td>%s</td>', all[names(columns)]), collapse='\n'),
        '</tr>'), collapse='\n')
}

#' Generate an formatted HTML table of parameter estimates
#' 
#' @param parframe A `data.frame` such as returned by [pmxparframe].
#' @param columns A named `character` vector of columns to include in the table
#' (and in what order). The names correspond to column names in `parframe` and
#' the value to the column labels that appear in the formatted table.
#' @param sections A `logical` indicating whether or not the table should be
#' formatted into sections according the the `type` column of `parframe`.
#' @param section.labels A named `character` vector. The names correspond to
#' values in the `type` column of `parframe`, and the values to labels that
#' appear in the formatted table.
#' @param footnote A `character` vector of footnotes to place underneath the
#' formatted table (may contain HTML codes).
#' @param show.fixed.to.zero A `logical` indicating whether parameters that are
#' fixed to zero should appear in the formatted table (by default, parameters
#' that are formatted to values other than zero do appear in the table, but
#' those that are fixed to zero are ignored).
#' @param merge.units A `logical` indicating whether or not units (if present)
#' should be merged into the parameter label (i.e., in parentheses following
#' the name/label).
#' @param na A `character` string to use in the formatted table to indicate
#' missing or non-applicable values.
#' @param digits Number of significant digits to include in the formatted
#' table.
#' 
#' @return An object of class `"pmxpartab"`. This is essentially just an HTML
#' character string that displays in the default web browser or viewer when
#' printed (as per [htmltools::print.html()]).
#' 
#' @seealso [pmxparframe]
#'
#' @examples
#' \donttest{
#' outputs <- list(
#'   est = list(
#'     th = list(CL = 0.482334, VC = 0.0592686),
#'     om = list(nCL = 0.315414, nVC = 0.536025),
#'     sg = list(ERRP = 0.0508497)),
#'   se = list(
#'     th = list(CL = 0.0138646, VC = 0.00555121),
#'     om = list(nCL = 0.0188891, nVC = 0.0900352),
#'     sg = list(ERRP = 0.00182851)),
#'   fixed = list(
#'     th = list(CL = FALSE, VC = FALSE),
#'     om = list(nCL = FALSE, nVC = FALSE),
#'     sg = list(ERRP = FALSE)),
#'   shrinkage = list(nCL = 9.54556, nVC = 47.8771))
#' 
#' meta <- list(
#'   parameters = list(
#'     list(name="CL", label="Clearance", units="L/h", type="Structural"),
#'     list(name="VC", label="Volume", units="L", type="Structural", trans="exp"),
#'     list(name="nCL", label="On Clearance", type="IIV", trans="SD (CV%)"),
#'     list(name="nVC", label="On Volume", type="IIV"),
#'     list(name="ERRP", label="Proportional Error", units="%", type="RUV", trans="%")))
#' 
#' pmxpartab(pmxparframe(outputs, meta),
#'     columns=c(est="Estimate", rse="RSE%", ci95="95% CI", shrinkage="Shrinkage"),
#'     footnote="CI=confidence interval; RSE=relative standard error.")
#' 
#' 
#' # An example using a Cox model, where we construct the parframe manually:
#' library(survival)
#' cph.fit <- coxph(Surv(time, status) ~ ph.ecog + age, data=lung)
#' parframe <- with(summary(cph.fit), data.frame(
#'     name  = c("ph.ecog", "age"),
#'     label = c("ECOG performance score", "Age"),
#'     est   = coefficients[,"exp(coef)"],
#'     pval  = coefficients[,"Pr(>|z|)"],
#'     lci95 = conf.int[,"lower .95"],
#'     uci95 = conf.int[,"upper .95"]
#' ))
#' pmxpartab(parframe=parframe,
#'     columns=c(est="HR", ci95="95% CI", pval="P-Value"))
#' }
#' @export
pmxpartab <- function(
    parframe,

    columns=c(est="Estimate", rse="RSE%", ci95="95% CI", shrinkage="Shrinkage"),

    sections = !is.null(parframe$type),
    section.labels = c(
        Structural      = "Typical Values",
        CovariateEffect = "Covariate Effects",
        RUV             = "Residual Error",
        IIV             = "Between Subject Variability",
        IOV             = "Inter-Occasion Variability"),

    footnote=NULL,

    show.fixed.to.zero = FALSE,
    merge.units        = TRUE,
    na                 = "-",
    digits             = 3) {

    if (isFALSE(show.fixed.to.zero) & !is.null(parframe$fixed)) {
        parframe <- parframe[!(parframe$fixed & parframe$est==0),, drop=FALSE]
    }

    if (is.null(parframe$shrinkage)) {
        columns <- columns[names(columns) != "shrinkage"]
    }

    ncolumns <- length(columns) + 1

    thead <- paste0('<tr>\n<th>Parameter</th>\n',
        paste0(paste0('<th>', columns, '</th>'), collapse="\n"), '\n</tr>\n')

    thead <- paste0('<thead>\n', thead, '\n</thead>\n')


    tbody <- ""
    for (i in 1:nrow(parframe)) {
        if (isTRUE(sections)) {
            newsection <- (!is.null(parframe$type) && !is.na(parframe$type[i]) && (i == 1 || parframe$type[i] != parframe$type[i-1]))
            if (isTRUE(newsection)) {
                type <- parframe$type[i]
                if (type %in% names(section.labels)) {
                    label <- section.labels[[type]]
                } else {
                    label <- type
                }

                tbody <- paste0(tbody, partab_section(label, ncolumns=ncolumns), '\n')
            }
        }
        args <- c(parframe[i,, drop=FALSE], list(merge.units=merge.units, na=na, digits=digits, indent=sections, columns=columns))
        tbody <- paste0(tbody, do.call(partab_row, args), '\n')
    }
    tbody <- paste0('<tbody>\n', tbody, '\n</tbody>\n')

    if (!is.null(footnote)) {
        footnote <- sprintf('<p>%s</p>\n', footnote)
        footnote <- paste0(footnote, collapse="\n")
        tfoot <- sprintf('<tfoot><tr><td colspan="%d" class="partabfootnote">%s</td></tr></tfoot>\n', ncolumns, footnote)
    } else {
        tfoot <- ""
    }

    table <- paste0('<table>\n', thead, tbody, tfoot, '</table>\n')
    structure(table, class=c("pmxpartab", "html", "character"), html=TRUE)
}


#' Print `pmxpartab` object
#'
#' @param x An object returned by [pmxpartab].
#' @param ... Further arguments passed on to other `print` methods.
#'
#' @return Returns `x` invisibly.
#'
#' @details In an interactive context, the rendered table will be displayed in
#' a web browser. Otherwise, the HTML code will be printed as text.
#' @export
print.pmxpartab <- function(x, ...) {
    if (interactive()) {
        z <- htmltools::HTML(x)
        default.style <- htmltools::htmlDependency("pmxpartab", "1.0",
            src=system.file(package="pmxpartab", "pmxpartab_defaults_1.0"),
            stylesheet="pmxpartab_defaults.css")
        z <- htmltools::div(class="Rpmxpartab", default.style, z)
        z <- htmltools::browsable(z)
        print(z, ...) # Calls htmltools:::print.html(z, ...)
    } else {
        cat(x)
    }
    invisible(x)
}


#' Method for printing in a `knitr` context
#'
#' @param x An object returned by [pmxpartab].
#' @param ... Further arguments passed on to [knitr::knit_print].
#'
#' @return A 'character` vector (see [knitr::knit_print]).
#'
#' @importFrom knitr knit_print
#' @export
knit_print.pmxpartab <- function(x, ...) {
    knit_to_html <-
        !is.null(knitr::opts_knit$get("rmarkdown.pandoc.to")) &&
        grepl("^html", knitr::opts_knit$get("rmarkdown.pandoc.to"))

    if (knit_to_html) {
        z <- htmltools::HTML(x)
        default.style <- htmltools::htmlDependency("pmxpartab", "1.0",
            src=system.file(package="pmxpartab", "pmxpartab_defaults_1.0"),
            stylesheet="pmxpartab_defaults.css")
        z <- htmltools::div(class="Rpmxpartab", default.style, z)
        knitr::knit_print(z, ...)
    } else {
        knitr::knit_print(as.character(x), ...)
    }
}

#' Parse a parameter description string
#'
#' @param string A `character`, starting with a name, followed by an optional
#' comma separated list of key-value pairs that can be parsed as R code. See
#' Examples.
#'
#' @return A named `list`.
#'
#' @examples
#' # Example 1: all elements present
#' x <- "CL, label='Clearance', units='L/h', trans=exp, type='Structural'"
#' parse_parameter_description(x)
#' 
#' # Example 2: Some elements missing (trans), will take default value (NULL)
#' x <- "CL, label='Clearance', units='L/h', type='Structural'"
#' parse_parameter_description(x)
#' 
#' # Example 3: Only the name is given
#' x <- "CL"
#' parse_parameter_description(x)
#' 
#' # Example 4: positional arguments
#' x <- "CL, 'Clearance', 'L/h', type='Structural'"
#' parse_parameter_description(x)
#' @export
parse_parameter_description <- function(string) {
    # Returns a structured object representing a description of a parameter
    # (in this example just a list with some attributes; only the name is mandatory)
    parameter_description <- function(name, label=NULL, units=NULL, trans=NULL, type=NULL, ...) {
        list(name=name, label=label, units=units, trans=trans, type=type, ...)
    }

    x <- str2lang(paste0("parameter_description(", string, ")"))
    x[[2]] <- as.character(x[[2]]) # Interpret the first element (name) as a string even if not quoted
    eval(x)
}


#' Format p-values
#'
#' @param pval A numeric vector of p-values.
#' @param digits The number of significant digits to retain.
#' @param eps A numeric value. Under this threshold, rather than showing the
#' p-value itself, show "< 1e-X" where X is the largest integer satisfying
#' this relationship.
#' @param alpha The significance level.
#' @param star.symbol A character to display next to those p-values that are
#' statistically significant (i.e., less then `alpha`).
#' @param html A logical flag indicating whether to return HTML code or plain
#' text.
#' @param unicode.le A logical flag indicating whether to use unicode
#' symbol [U+2264](https://www.compart.com/en/unicode/U+2264)
#' for "less-than-or-equal-to" (only applies when `html` is `FALSE`).
#' @return A character vector of the same length as `pval`.
#' @seealso [base::format.pval]
#' @examples
#' x <- c(1, 0.5, 0.05, 0.049, 0.01, 0.001, 0.0001, 0.00001)
#' fpval(x, html=FALSE, unicode.le=FALSE)
#' @export
fpval <- function(pval, digits=3, eps=1e-3, alpha=0.05, star.symbol="*", html=FALSE, unicode.le=FALSE) {
    .fpval.internal <- function(pval) {
        if (pval >= eps) {
            pval <- format.pval(pval, digits=digits, eps=eps)
            if (substring(pval, 1, 1) == "<") {
                if (html) {
                    pval <- paste0("&lt; ", substring(pval, 2))
                } else {
                    pval <- paste0("< ", substring(pval, 2))
                }
            } else {
                pval <- paste0("", pval)
            }
        } else {
            if (html) {
                pval <- paste0("&le; 10<sup>", ceiling(log10(pval)), "</sup>")
            } else if (unicode.le) {
                pval <- paste0("\U{2264} 1e", ceiling(log10(pval)))
            } else {
                pval <- paste0("<= 1e", ceiling(log10(pval)))
            }
        }
        pval
    }
    star <- ifelse(pval < alpha, star.symbol, "")
    paste0(Vectorize(.fpval.internal)(pval), star)
}

