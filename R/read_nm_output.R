#' Read meta information from a YAML file
#'
#' @param meta.file The name of a YAML file to be read.
#' @return A `list` if the file exists, otherwise `NULL`.
#' @export
read_meta <- function(meta.file="meta.yaml") {
    if (file.exists(meta.file)) {
        yaml::read_yaml(meta.file)
    } else {
        NULL
    }
}

#' Read NONMEM output
#'
#' @details All arguments are optional. If a particular output file cannot be
#' found, then it is simply skipped (and the resulting object won't contain the
#' components that would normally be read from there).
#'
#' @param rundir Name of the directory containing the output files.
#' @param runname Name of the run (i.e., corresponds to the basename of the output files).
#' @param lst.file Name of the .lst file (standard NONMEM output file).
#' @param ext.file Name of the .ext file (standard NONMEM output file).
#' @param shk.file Name of the .shk file (standard NONMEM output file).
#' @param phi.file Name of the .phi file (standard NONMEM output file).
#' @param phm.file Name of the .phm file (standard NONMEM output file).
#' @param cov.file Name of the .cov file (standard NONMEM output file).
#' @param cor.file Name of the .cor file (standard NONMEM output file).
#' @param bootstrap.file Name of the file containing bootstrap results (typically produced by PsN).
#' @param meta Object containing meta information that accompanies the model
#' (e.g., names of the `THETA`, `OMEGA` and `SIGMA` parameters).
#' @param th_names A character vector containing the names associated with the
#' `THETA` parameters in their respective order (e.g., `THETA(1)` is given the
#' name `th_names[1]`, and so on).
#' @param om_names A character vector containing the names associated with the
#' the `OMEGA` matrix diagonal elements, in their respective order (e.g.,
#' `OMEGA(1,1)` is given the name `om_names[1]`, and so on).
#' @param sg_names A character vector containing the names associated with the
#' `SIGMA` matrix diagonal elements, in their respective order (e.g.,
#' `SIGMA(1,1)` is given the name `sg_names[1]`, and so on).
#' @param use.vcov Should the default `OMEGA` and `SIGMA` be on the
#' variance/covariance scale instead of the SD/correlation scale?
#' @param ... Additional arguments (ignored).
#'
#' @return A named list with components containing the outputs from the NONMEM
#' run. Notably, the components `th`, `om` and `sg` contain the final estimates
#' of the `THETA`, `SD(ETA)` and `SD(EPS)` parameters respectively (`SD` means
#' standard deviation).
#'
#' @importFrom utils read.csv read.table tail
#' @importFrom stats median quantile sd setNames
#' @export
read_nm_output <- function(
    rundir       = getwd(),
    runname      = basename(normalizePath(rundir)),
    lst.file     = file.path(rundir, sprintf("%s.lst", runname)),
    ext.file     = file.path(rundir, sprintf("%s.ext", runname)),
    shk.file     = file.path(rundir, sprintf("%s.shk", runname)),
    phi.file     = file.path(rundir, sprintf("%s.phi", runname)),
    phm.file     = file.path(rundir, sprintf("%s.phm", runname)),
    cov.file     = file.path(rundir, sprintf("%s.cov", runname)),
    cor.file     = file.path(rundir, sprintf("%s.cor", runname)),
    bootstrap.file = file.path(rundir, "bootstrap", sprintf("raw_results_%s.csv", runname)),
    meta         = read_meta(file.path(rundir, "meta.yaml")),
    th_names     = meta$namemap$theta,
    om_names     = meta$namemap$omega,
    sg_names     = meta$namemap$sigma,
    use.vcov     = FALSE,
    ...) {

    res <- list()
    res$meta <- meta

    # Read .ext file
    if (is.character(ext.file) && file.exists(ext.file)) {
        l <- readLines(ext.file)
        x <- grep("TABLE NO\\.", l)
        ext <- read.table(text=paste(l[ (x[length(x)]) : length(l) ], collapse="\n"), skip=1, header=TRUE)
        names(ext) <- tolower(names(ext))

        LTmat <- function (LT) 
        {
            x <- length(LT)
            p <- (sqrt(8 * x + 1) - 1)/2
            m <- matrix(0, p, p)
            m[upper.tri(m, diag = T)] <- LT
            m2 <- t(m)
            diag(m2) <- 0
            m + m2
        }

        if (names(ext)[ncol(ext)] != "obj") {
            warning("This form of estimation is not supported at the moment.")
        }

        ifinal <- NULL
        ise <- NULL
        icor <- NULL
        isecor <- NULL
        ifixed <- NULL
        if (any(ext$iteration==-1000000000)) {
            ifinal <- tail(which(ext$iteration==-1000000000), 1)
        }
        if (any(ext$iteration==-1000000001)) {
            ise <- tail(which(ext$iteration==-1000000001), 1)
        }
        if (any(ext$iteration==-1000000004)) {
            icor <- tail(which(ext$iteration==-1000000004), 1)
        }
        if (any(ext$iteration==-1000000005)) {
            isecor <- tail(which(ext$iteration==-1000000005), 1)
        }
        if (any(ext$iteration==-1000000006)) {
            ifixed <- tail(which(ext$iteration==-1000000006), 1)
        }
        if (!is.null(ifinal)) {

            ofv <- ext[ifinal, grepl("obj", names(ext))]
            ofv <- as.numeric(ofv)

            th <- ext[ifinal, grepl("^theta", names(ext))]
            th <- as.numeric(th)
            names(th) <- th_names

            om_cov <- ext[ifinal, grepl("^omega", names(ext))]
            om_cov <- LTmat(as.numeric(om_cov))
            dimnames(om_cov) <- list(om_names, om_names)

            sg_cov <- ext[ifinal, grepl("^sigma", names(ext))]
            sg_cov <- LTmat(as.numeric(sg_cov))
            dimnames(sg_cov) <- list(sg_names, sg_names)

            if (!is.null(icor)) {

                om_cor <- ext[icor, grepl("^omega", names(ext))]
                om_cor <- LTmat(as.numeric(om_cor))
                dimnames(om_cor) <- list(om_names, om_names)

                sg_cor <- ext[icor, grepl("^sigma", names(ext))]
                sg_cor <- LTmat(as.numeric(sg_cor))
                dimnames(sg_cor) <- list(sg_names, sg_names)
            }

            if (use.vcov) {
                om <- diag(om_cov)
                sg <- diag(sg_cov)
            } else {
                om <- diag(om_cor)
                sg <- diag(sg_cor)
            }

            all <- c(th, om, sg)

            res$ofv    <- ofv
            res$all    <- all
            res$th     <- th
            res$om     <- om
            res$sg     <- sg
            res$om_cov <- om_cov
            res$sg_cov <- sg_cov
            res$om_cor <- om_cor
            res$sg_cor <- sg_cor
        }

        if (!is.null(ifixed)) {
            th_fix  <- ext[ifixed, grepl("^theta", names(ext))]
            th_fix <- as.numeric(th_fix) == 1
            names(th_fix) <- th_names

            om_cov_fix <- ext[ifixed, grepl("^omega", names(ext))]
            om_cov_fix <- LTmat(as.numeric(om_cov_fix)) == 1
            dimnames(om_cov_fix) <- list(om_names, om_names)

            sg_cov_fix <- ext[ifixed, grepl("^sigma", names(ext))]
            sg_cov_fix <- LTmat(as.numeric(sg_cov_fix)) ==1
            dimnames(sg_cov_fix) <- list(sg_names, sg_names)

            om_cor_fix <- ext[ifixed, grepl("^omega", names(ext))]
            om_cor_fix <- LTmat(as.numeric(om_cor_fix)) == 1
            dimnames(om_cor_fix) <- list(om_names, om_names)

            sg_cor_fix <- ext[ifixed, grepl("^sigma", names(ext))]
            sg_cor_fix <- LTmat(as.numeric(sg_cor_fix)) ==1
            dimnames(sg_cor_fix) <- list(sg_names, sg_names)

            if (use.vcov) {
                om_fix <- diag(om_cov_fix)
                sg_fix <- diag(sg_cov_fix)
            } else {
                om_fix <- diag(om_cor_fix)
                sg_fix <- diag(sg_cor_fix)
            }

            all_fix <- c(th_fix, om_fix, sg_fix)

            res$fixed$all    <- all_fix
            res$fixed$th     <- th_fix
            res$fixed$om     <- om_fix
            res$fixed$sg     <- sg_fix
            res$fixed$om_cov <- om_cov_fix
            res$fixed$sg_cov <- sg_cov_fix
            res$fixed$om_cor <- om_cor_fix
            res$fixed$sg_cor <- sg_cor_fix
        }

        if (!is.null(ise)) {
            th_se  <- ext[ise, grepl("^theta", names(ext))]
            th_se <- as.numeric(th_se)
            th_se[th_fix] <- NA
            names(th_se) <- th_names

            om_cov_se <- ext[ise, grepl("^omega", names(ext))]
            om_cov_se <- LTmat(as.numeric(om_cov_se))
            om_cov_se[om_cov_fix] <- NA
            dimnames(om_cov_se) <- list(om_names, om_names)

            sg_cov_se <- ext[ise, grepl("^sigma", names(ext))]
            sg_cov_se <- LTmat(as.numeric(sg_cov_se))
            sg_cov_se[sg_cov_fix] <- NA
            dimnames(sg_cov_se) <- list(sg_names, sg_names)

            if (!is.null(isecor)) {
                om_cor_se <- ext[isecor, grepl("^omega", names(ext))]
                om_cor_se <- LTmat(as.numeric(om_cor_se))
                om_cor_se[om_fix] <- NA
                dimnames(om_cor_se) <- list(om_names, om_names)

                sg_cor_se <- ext[isecor, grepl("^sigma", names(ext))]
                sg_cor_se <- LTmat(as.numeric(sg_cor_se))
                sg_cor_se[sg_fix] <- NA
                dimnames(sg_cor_se) <- list(sg_names, sg_names)
            }

            if (use.vcov) {
                om_se <- diag(om_cov_se)
                sg_se <- diag(sg_cov_se)
            } else {
                om_se <- diag(om_cor_se)
                sg_se <- diag(sg_cor_se)
            }

            all_se <- c(th_se, om_se, sg_se)

            res$se$all    <- all_se
            res$se$th     <- th_se
            res$se$om     <- om_se
            res$se$sg     <- sg_se
            res$se$om_cov <- om_cov_se
            res$se$sg_cov <- sg_cov_se
            res$se$om_cor <- om_cor_se
            res$se$sg_cor <- sg_cor_se
        }

    }

    # Read .shk file
    if (is.character(shk.file) && file.exists(shk.file)) {
        l <- readLines(shk.file)
        x <- grep("TABLE NO\\.", l)
        shk <- read.table(text=paste(l[ (x[length(x)]) : length(l) ], collapse="\n"), skip=1, header=TRUE)
        names(shk) <- tolower(names(shk))

        if (length(unique(shk$subpop)) == 1) {
            shrinkage <- shk[shk$type==4, grepl("^eta", names(shk))]
            shrinkage <- as.numeric(shrinkage)
            names(shrinkage) <- om_names
            res$shrinkage <- shrinkage
        } else {
            shrinkage <- lapply(split(shk, shk$subpop), function(.) {
                shrinkage <- .[.$type==4, grepl("^eta", names(.))]
                setNames(as.numeric(shrinkage), om_names)
            })
            res$mixture$shrinkage <- shrinkage
        }
    }

    # Read .lst file
    if (is.character(lst.file) && file.exists(lst.file)) {
        l <- readLines(lst.file)

        runstarted <- strptime(paste(l[1], l[2]), "%A %m/%d/%Y %I:%M %p")

        i <- which(grepl("^NM-TRAN MESSAGES", l))
        l <- l[(1:length(l))>i]

        convergence <- "FAILED"
        if (any(grepl("MINIMIZATION SUCCESSFUL", l))) {
            convergence <- "SUCCESSFUL"
        }
        if (any(grepl("HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION", l))) {
            convergence <- "PROBLEMS"
        }
        if (any(grepl("OPTIMIZATION WAS COMPLETED", l))) {
            convergence <- "SUCCESSFUL"
        }
        covstep <- convergence > 0 && !is.null(res$se)

        nmversion <- NULL
        i <- grepl("1NONLINEAR MIXED EFFECTS MODEL PROGRAM \\(NONMEM\\) VERSION", l)
        if (any(i)) {
            nmversion <- gsub("1NONLINEAR MIXED EFFECTS MODEL PROGRAM \\(NONMEM\\) VERSION", "", l[i])
            nmversion <- gsub("^\\s+", "", nmversion)
            nmversion <- gsub("\\s+$", "", nmversion)
        }

        n.obs <- NULL
        i <- grepl("TOT. NO. OF OBS RECS:", l)
        if (any(i)) {
            n.obs <- scan(quiet=T, text=gsub(".*:", "", l[i]))
        }

        n.indiv <- NULL
        i <- grepl("TOT. NO. OF INDIVIDUALS:", l)
        if (any(i)) {
            n.indiv <- scan(quiet=T, text=gsub(".*:", "", l[i]))
        }

        runtime.estim <- NULL
        i <- grepl("Elapsed estimation  time in seconds:", l)
        if (any(i)) {
            runtime.estim <- scan(quiet=T, text=gsub(".*:", "", l[i]))
        }

        runtime.covstep <- NULL
        i <- grepl("Elapsed covariance  time in seconds:", l)
        if (any(i)) {
            runtime.covstep <- scan(quiet=T, text=gsub(".*:", "", l[i]))
        }

        runtime.postproc <- NULL
        i <- grepl("Elapsed postprocess time in seconds:", l)
        if (any(i)) {
            runtime.postproc <- scan(quiet=T, text=gsub(".*:", "", l[i]))
        }

        res$nmversion        <- nmversion
        res$runstarted       <- runstarted
        res$runtime$estim    <- runtime.estim
        res$runtime$covstep  <- runtime.covstep
        res$runtime$postproc <- runtime.postproc
        res$num$indiv        <- n.indiv
        res$num$obs          <- n.obs
        res$minimization     <- convergence
        res$covstep          <- covstep
    }

    # Read .cov file
    if (is.character(cov.file) && file.exists(cov.file)) {
        l <- readLines(cov.file)
        x <- grep("TABLE NO\\.", l)
        cov <- read.table(text=paste(l[ (x[length(x)]) : length(l) ], collapse="\n"), skip=1, header=TRUE)
        names(cov) <- tolower(names(cov))
        cov <- cov[, -1]

        #if (!is.null(th_names)) {
        #    names(cov)[grepl("^theta", names(cov))] <- th_names
        #}
        #if (!is.null(sg_names)) {
        #    names(cov)[grepl("^sigma", names(cov))] <- sg_names
        #}
        #if (!is.null(om_names)) {
        #    names(cov)[grepl("^omega", names(cov))] <- om_names
        #}

        cov <- as.matrix(cov)
        rownames(cov) <- colnames(cov)
        res$se_cov <- cov
    }

    # Read .cor file
    if (is.character(cor.file) && file.exists(cor.file)) {
        l <- readLines(cor.file)
        x <- grep("TABLE NO\\.", l)
        cor <- read.table(text=paste(l[ (x[length(x)]) : length(l) ], collapse="\n"), skip=1, header=TRUE)
        names(cor) <- tolower(names(cor))
        cor <- cor[, -1]

        #if (!is.null(th_names)) {
        #    names(cor)[grepl("^theta", names(cor))] <- th_names
        #}
        #if (!is.null(sg_names)) {
        #    names(cor)[grepl("^sigma", names(cor))] <- sg_names
        #}
        #if (!is.null(om_names)) {
        #    names(cor)[grepl("^omega", names(cor))] <- om_names
        #}

        cor <- as.matrix(cor)
        rownames(cor) <- colnames(cor)
        res$se_cor <- cor

        # Eigenvalues and condition number
        temp <- cor
        i <- apply(temp, 1, sd) > 0
        temp <- temp[i, i]
        diag(temp) <- 1
        eigv <- rev(eigen(temp)$values)  # From smallest to largest, the way NONMEM shows them in .lst
        res$eigv <- eigv
        res$condition_number <- max(eigv)/min(eigv)
    }

    # Read .phi file
    if (is.character(phi.file) && file.exists(phi.file)) {
        l <- readLines(phi.file)
        x <- grep("TABLE NO\\.", l)
        phi <- read.table(text=paste(l[ (x[length(x)]) : length(l) ], collapse="\n"), skip=1, header=TRUE)
        names(phi) <- tolower(names(phi))
        i <- grepl("^eta", names(phi))
        if (sum(i) == length(om_names)) {
            res$etas <- setNames(phi[,i], om_names)
        }
    }

    # Read .phm file (only for mixture models)
    if (is.character(phm.file) && file.exists(phm.file)) {
        l <- readLines(phm.file)
        x <- grep("TABLE NO\\.", l)
        phm <- read.table(text=paste(l[ (x[length(x)]) : length(l) ], collapse="\n"), skip=1, header=TRUE)
        names(phm) <- tolower(names(phm))
        i <- grepl("^eta", names(phm))
        if (sum(i) == length(om_names)) {
            names(phm)[i] <- om_names
        }

        if (length(unique(phm$subpop)) > 1) {
            res$mixture$nsubpop <- length(unique(phm$subpop))
            res$mixture$subpop <- split(phm, phm$subpop)
        }
    }

    # Read bootstrap results
    if (is.character(bootstrap.file) && file.exists(bootstrap.file)) {
        bootstrap.data <- read.csv(bootstrap.file, header=T, check.names=F)

        bs.ofv <- function(x) {
            is.ofv <- toupper(names(x)) == "OFV"
            x[, is.ofv, drop=F]
        }

        bs.th <- function(x) {
            is.th <- names(x) %in% th_names | grepl("^THETA", names(x))
            x[, is.th, drop=F]
        }

        bs.matrix <- function(x, matrx=c("OMEGA", "SIGMA")) {

            nx <- names(x)
            names(nx) <- nx

            if (matrx == "OMEGA") {
                d <- length(res$om)
                q <- om_names
            } else {
                d <- length(res$sg)
                q <- sg_names
            }
            nm <- outer(1:d, 1:d, function(x, y) paste0(matrx, "(", y, ",", x, ")"))
            dnm <- diag(nm)
            nx[q] <- dnm
            nm <- nm[upper.tri(nm, diag=T)]
            names(x) <- nx
            m <- matrix(0, nrow=nrow(x), ncol=length(nm))
            m <- as.data.frame(m)
            names(m) <- nm

            i <- intersect(nm, nx)
            m[,i] <- x[,i]

            if (!use.vcov) {
                # variance-covariance to sd/cor
                if (ncol(m) == 1) {
                    m <- sqrt(m)
                } else {
                    m <- t(apply(m, 1, function(x) {
                        #if (length(x) == 1) return(sqrt(x))
                        v <- LTmat(x)
                        s <- diag(1/sqrt(diag(v)))
                        r <- s %*% v %*% s
                        diag(r) <- sqrt(diag(v))
                        r[is.nan(r)] <- 0
                        t(r)[upper.tri(r, diag=T)]
                    }))
                    m <- as.data.frame(m)
                }
            }
            names(nm) <- nm
            nm[dnm] <- q
            names(m) <- nm
            m
        }

        bs.om <- function(x) { bs.matrix(x, "OMEGA") }
        bs.sg <- function(x) { bs.matrix(x, "SIGMA") }

        bootstrap.keep <- cbind(
            bs.ofv(bootstrap.data),
            bs.th(bootstrap.data),
            bs.om(bootstrap.data),
            bs.sg(bootstrap.data))

        boot.fixed <- sapply(bootstrap.keep, function(x) length(unique(x[!is.na(x)]))) == 1
        bootstrap.keep <- bootstrap.keep[, !boot.fixed, drop=F]

        bootstrap.orig <- bootstrap.keep[1,, drop=F]
        bootstrap.keep <- bootstrap.keep[-1,, drop=F]
        bootstrap.data <- bootstrap.data[-1,, drop=F]

        success <- bootstrap.data$minimization_successful == 1

        res$bootstrap$convergence <- bootstrap.data[,
            c("minimization_successful", "covariance_step_successful",
                "estimate_near_boundary", "rounding_errors")]
        res$bootstrap$n$total          <- nrow(bootstrap.data)
        res$bootstrap$n$successful     <- sum(bootstrap.data$minimization_successful, na.rm=T)
        res$bootstrap$n$covstep        <- sum(bootstrap.data$covariance_step_successful, na.rm=T)
        res$bootstrap$n$nearboundary   <- sum(bootstrap.data$estimate_near_boundary, na.rm=T)
        res$bootstrap$n$roundingerrors <- sum(bootstrap.data$rounding_errors, na.rm=T)

        res$bootstrap$data        <- bootstrap.keep
        res$bootstrap$orig        <- bootstrap.orig
        res$bootstrap$median      <- sapply(res$bootstrap$data[success,], median, na.rm=T)
        res$bootstrap$ci          <- lapply(res$bootstrap$data[success,], quantile, probs=c(0.025, 0.975), na.rm=T)
        res$bootstrap$ci          <- as.data.frame(res$bootstrap$ci, optional=T)
        res$bootstrap$bias        <- 100*(res$bootstrap$median - res$bootstrap$orig)/res$bootstrap$orig

    }

    res
}


