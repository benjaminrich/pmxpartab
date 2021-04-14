#' Get parameters in a data.frame
#'
#' @export
parframe <- function(outputs, meta=outputs$meta) {
    param <- meta$parameters
    z <- data.table::rbindlist(param, fill=T)

    outputs$all       <- with(outputs, c(th, om, sg))
    outputs$se$all    <- with(outputs$se, c(th, om, sg))
    outputs$fixed$all <- with(outputs$fixed, c(th, om, sg))

    z$fixed     <- as.logical(NA)
    z$est       <- as.numeric(NA)
    z$se        <- as.numeric(NA)
    z$rse       <- as.numeric(NA)
    z$lci95     <- as.numeric(NA)
    z$uci95     <- as.numeric(NA)
    z$pval      <- as.numeric(NA)
    z$shrinkage <- as.numeric(NA)

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
        if (name %in% names(outputs$all)) {
            est <- outputs$all[[name]]
            se <- outputs$se$all[[name]]
            fixed <- outputs$fixed$all[[name]]
            if (have.bootstrap && name %in% names(outputs$bootstrap$median)) {
                boot.median <- outputs$bootstrap$median[[name]]
                boot.lci <- outputs$bootstrap$ci95[[name]][1]
                boot.uci <- outputs$bootstrap$ci95[[name]][2]
            }
        } else if ((re <- regexec("^om(?:_cor)?\\((\\w+),(\\w+)\\)$", name, perl=T))[[1]][1] > 0) {
            m1 <- regmatches(name, re)[[1]][2]
            m2 <- regmatches(name, re)[[1]][3]
            om_names <- names(outputs$om)
            a1 <- match(m1, om_names)
            a2 <- match(m2, om_names)
            if (!is.na(a1) && !is.na(a2)) {
                est <- outputs$om_cor[a1,a2]
                se <- outputs$se$om_cor[a1,a2]
                fixed <- outputs$fixed$om_cor[a1,a2]
                if (fixed) {
                    se <- NA
                }
                if (have.bootstrap) {
                    boot.re <- paste0("OMEGA\\D", max(a1, a2), "\\D", min(a1, a2), "($|\\D)")
                    boot.name <- grep(boot.re, names(nmout$bootstrap$median), ignore.case=T, value=T)
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
            if (a1 %in% dimnames(outputs$om_cov)[[1]] && a2 %in% dimnames(outputs$om_cov)[[2]]) {
                est <- outputs$om_cov[a1,a2]
                se <- outputs$se$om_cov[a1,a2]
                fixed <- outputs$fixed$om_cov[a1,a2]
                if (fixed) {
                    se <- NA
                }
                if (have.bootstrap) {
                    boot.re <- paste0("OMEGA\\D", max(a1, a2), "\\D", min(a1, a2), "($|\\D)")
                    boot.name <- grep(boot.re, names(nmout$bootstrap$median), ignore.case=T, value=T)
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

        if (fixed || is.null(se) || is.na(se)) {
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

        if (z$name[i] %in% names(outputs$shrinkage)) {
            z$shrinkage[i] <- outputs$shrinkage[z$name[i]]
        }
    }
    z <- subset(z, !is.na(est))
    as.data.frame(z)
}

# Internal function to help format numbers
p <- function(x, digits=3, flag="", round.integers=FALSE){
    if (!is.numeric(x)) {
        return(x)
    }
    prefix <- ifelse(flag=="+" & x > 0, "+", "")
    paste0(prefix, table1::signif_pad(x, digits=digits, round.integers=round.integers))
}

parameter.estimate.table.section <- function(label, ncolumns) {
    paste0(c('<tr>',
        paste0(sprintf('<td class="paramsectionheading">%s</td>', c(label, rep("", ncolumns-1))), collapse='\n'),
        '</tr>'), collapse='\n')
}

parameter.estimate.table.row <- function(
    name,
    label          = NULL,
    units          = NULL,
    type           = c("Structural", "CovariateEffect", "IIV", "IOV", "RUV", "Unspecified"),
    trans          = c("identity", "%", "exp", "ilogit", "CV%", "SD (CV%)"),
    expression     = NULL,
    relatedTo      = NULL,
    superscript    = NULL,
    fixed          = NULL,
    est            = NULL,
    se             = NULL,
    rse            = NULL,
    lci95          = NULL,
    uci95          = NULL,
    boot.median    = NULL,
    boot.lci       = NULL,
    boot.uci       = NULL,
    shrinkage      = NULL,
    na             = "n/a",
    digits         = 3,
    indent         = TRUE,
    have.bootstrap = !is.null(boot.median),
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
        label <- sprintf("%s (%s)", label, units)
    }

    if (!is.null(trans) && !is.na(trans) && trans == "SD (CV%)") {
        g <- function(x) { 100*sqrt(exp(x^2) - 1) }
        x <- est
        est <- sprintf("%s (%s%%)", p(x, digits), p(g(x), digits))
    } else {
        est <- p(est, digits)
    }

    est <- paste0(est, superscript)

    if (fixed) {
        est <- sprintf('%s Fixed', est)
    }

    if (is.na(se)) {
        se <- na
        rse <- na
        ci95 <- na
    } else {
        rse <- p(rse, digits)
        ci95 <- sprintf('%s &ndash; %s', p(lci95, digits), p(uci95, digits))
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
            shrinkage <- ""
        } else {
            shrinkage <- sprintf("%s%%", p(shrinkage, digits))
        }
    }
    all <- c(est=est, rse=rse, ci95=ci95, boot.median=boot.median, boot.ci95=boot.ci95, shrinkage=shrinkage)
    paste0(c('<tr>',
        sprintf('<td class="%s">%s</td>', ifelse(isTRUE(indent), "paramlabelindent", "paramlabelnoindent"), label),
        paste0(sprintf('<td>%s</td>', all[names(columns)]), collapse='\n'),
        '</tr>'), collapse='\n')
}

#' Generate a parameter estimates table in HTML
#' 
#' @examples
#' outputs <- list(
#'     th = list(CL = 0.482334, VC = 0.0592686),
#'     om = list(nCL = 0.315414, nVC = 0.536025),
#'     sg = list(ERRP = 0.0508497),
#'     se = list(
#'         th = list(CL = 0.0138646, VC = 0.00555121),
#'         om = list(nCL = 0.0188891, nVC = 0.0900352),
#'         sg = list(ERRP = 0.00182851)),
#'     fixed = list(
#'         th = list(CL = FALSE, VC = FALSE),
#'         om = list(nCL = FALSE, nVC = FALSE),
#'         sg = list(ERRP = FALSE)),
#'     shrinkage = list(nCL = 9.54556, nVC = 47.8771))
#' 
#' meta <- list(
#'     parameters = list(
#'         list(name="CL",   label="Clearance",    units="L/h", type="Structural"),
#'         list(name="VC",   label="Volume",       units="L",   trans="exp", type="Structural"),
#'         list(name="nCL",  label="On Clearance",              trans="SD (CV%)", type="IIV"),
#'         list(name="nVC",  label="On Volume",                 type="IIV"),
#'         list(name="ERRP", label="Proportional Error",        units="%", trans="%", type="RUV")))
#' 
#' parframe(outputs, meta)
#' 
#' pmxpartab(parframe(outputs, meta),
#'     columns=c(est="Estimate", rse="RSE%", ci95="95%CI", shrinkage="Shrinkage"))
#' @export
pmxpartab <- function(
    parframe,

    columns=c(est="Estimate", rse="RSE%", ci95="95%CI", shrinkage="Shrinkage"),

    sections = TRUE,
    section.labels = c(
        Structural      = "Typical Values",
        CovariateEffect = "Covariate Effects",
        RUV             = "Residual Error",
        IIV             = "Between Subject Variability",
        IOV             = "Inter-Occasion Variability"),

    show.fixed.to.zero=F,
    na="n/a",
    digits=3) {

    if (isFALSE(show.fixed.to.zero)) {
        parframe <- subset(parframe, !(fixed & est==0))
    }

    ncolumns <- length(columns) + 1

    thead <- paste0('<tr>\n<th rowspan="2">Parameter</th>\n',
        paste0(paste0('<th rowspan="2">', columns, '</th>'), collapse="\n"), '\n</tr>\n')


    tbody <- ""
    for (i in 1:nrow(parframe)) {
        if (isTRUE(sections)) {
            newsection <- (!is.null(parframe$type) && !is.na(parframe$type[i]) && (i == 1 || parframe$type[i] != parframe$type[i-1]))
            if (newsection) {
                type <- parframe$type[i]
                if (type %in% names(meta$labels)) {
                    label <- meta$labels[[type]]
                } else if (type %in% names(section.labels)) {
                    label <- section.labels[[type]]
                } else {
                    label <- type
                }

                tbody <- paste0(tbody, parameter.estimate.table.section(label, ncolumns=ncolumns), '\n')
            }
        }
        args <- c(parframe[i,], list(na=na, digits=digits, indent=sections))
        tbody <- paste0(tbody, do.call(parameter.estimate.table.row, args), '\n')
    }

    table <- paste0('<table>\n<thead>\n', thead, '\n</thead>\n<tbody>\n', tbody, '\n</tbody>\n</table>\n')
    structure(table, class=c("pmxpartab", "html", "character"), html=TRUE)
}


#' Print \code{pmxpartab} object.
#' @param x An object returned by \code{\link{pmxpartab}}.
#' @param ... Further arguments passed on to other \code{print} methods.
#' @return Returns \code{x} invisibly.
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


#' Method for printing in a \code{knitr} context.
#' @param x An object returned by \code{\link{pmxpartab}}.
#' @param ... Further arguments passed on to \code{knitr::knit_print}.
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
    parameter_description <- function(name, label=NULL, units=NULL, trans=NULL, type=NULL) {
        list(name=name, label=label, units=units, trans=trans, type=type)
    }

    x <- str2lang(paste0("parameter_description(", string, ")"))
    x[[2]] <- as.character(x[[2]]) # Interpret the first element (name) as a string even if not quoted
    eval(x)
}

