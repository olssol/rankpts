## -----------------------------------------------------------------------------
##
##  DESCRIPTION:
##      This file contains R functions for dose selection fix studies
##
##  DATE:
##      AUGUST, 2023
## -----------------------------------------------------------------------------

#' Describe the design
#'
#' @export
#'
desp1_describe <- function(x, ...) {
    cat("Type:\n")
    cat("    Phase I dose escaltion studies \n\n")
    cat("Design Parameters:\n")
    cat("    sample_size:   maximum sample size of the study \n")
    cat("                   (default 24)\n")
    cat("    size_dose:     maximum sample size of each dose level \n")
    cat("                   (default 9)\n")
    cat("    size_cohort:   cohort size \n")
    cat("                   (default 3)\n")
    cat("    tox_rate:      toxicity rates \n")
    cat("                   (default c(0.1, 0.2, 0.3))\n")
    cat("    res_rate:      response rates \n")
    cat("                   (default c(0.2, 0.4, 0.6))\n")
    cat("    rho:           odds ratio \n")
    cat("                   (default 1)\n")
    cat("    dlt_days:      number of days of the DLT window \n")
    cat("                   (default 28)\n")
    cat("    intra_fac:     factor for toxicity rate if intra dose
                            escalation \n")
    cat("                   (default 1.3) \n")
    cat("    par_enroll:    list of enrollment parameters \n")
}


#' Default design parameter
#'
#'
internal_desp1_dpara <- function() {
    rst <- list(sample_size       = 24,
                size_dose         = 9,
                size_cohort       = 3,
                tox_rate          = c(0.1, 0.2, 0.3),
                res_rate          = c(0.2, 0.4, 0.6),
                rho               = 1,
                dlt_days          = 28,
                intra_fac         = 1.3,
                par_enroll        = list(type       = "by_rate",
                                         pt_per_mth = 3)
                )

    rst$theta <- stb_tl_p1_sce(rst$tox_rate, rst$res_rate, rst$rho)
    rst
}

#' Generate data
#'
#'
#'
desp1_gen_data <- function(lst_design, seed = NULL, ...) {

    if (!is.null(seed))
        old_seed <- set.seed(seed)

    rst <- NULL

    ## reset
    if (!is.null(seed))
        set.seed(old_seed)

    ## return
    rst
}
