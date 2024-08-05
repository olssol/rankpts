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
    cat("Design Parameters:\n")
    cat("    sample_size:   maximum sample size of the study \n")
    cat("                   (default 24)\n")
    cat("    size_dose:     maximum sample size of each dose level \n")
    cat("                   (default 9)\n")
    cat("    tox_rate:      toxicity rates \n")
    cat("                   (default c(0.1, 0.2, 0.3))\n")
    cat("    res_rate:      response rates \n")
    cat("                   (default c(0.2, 0.4, 0.6))\n")
    cat("    rho:           odds ratio \n")
    cat("                   (default 1)\n")
    cat("    dlt_days:      number of days of the DLT window \n")
    cat("                   (default 28)\n")
    cat("    size_cohort_regular: regular cohort size without accelaration \n")
    cat("                   (default 3)\n")
    cat("    n_reuse_regular: number of titrated patients in \n")
    cat("                   regular cohorts  (default 0)\n")
    cat("    size_cohort_acc: acceleration phase cohort size\n")
    cat("                   (default 1)\n")
    cat("    n_reuse_acc:   number of titrated patients in acceleration \n")
    cat("                   (default 1)\n")
    cat("    acc_max_dose:  maximum number of doses for titrated patients \n")
    cat("                   (default 3)\n")
    cat("    intra_fac:     factor for toxicity rate if intra dose\n")
    cat("                   escalation (default 1.1) \n")
    cat("    par_enroll:    list of enrollment parameters \n")
}

#' Default design parameter
#'
inter_desp1_dpara_ext <- function(rst) {
    rst$theta  <- stb_tl_p1_sce(rst$tox_rate, rst$res_rate, rst$rho)
    rst$n_dose <- length(rst$tox_rate)
    rst
}

#' Default design parameter
#'
#'
internal_desp1_dpara <- function() {
    rst <- list(sample_size       = 24,
                size_dose         = 9,
                tox_rate          = c(0.1, 0.2, 0.3),
                res_rate          = c(0.2, 0.4, 0.6),
                rho               = 1,
                dlt_days          = 28,
                size_cohort_regular = 3,
                n_reuse_regular     = 0,
                size_cohort_acc     = 3,
                n_reuse_acc         = 0,
                acc_max_dose        = 3,
                intra_fac           = 1.1,
                par_enroll          = list(type       = "by_rate",
                                           pt_per_mth = 1),
                date_bos          = "2024-01-01"
                )

    inter_desp1_dpara_ext(rst)
}

#' Generate Cohort
#'
#'
#'
desp1_generate_cohort_flex <- function(lst_para, dose, data, ...) {

    ## regular stage
    size_cohort_regular <- lst_para$size_cohort_regular
    n_reuse_regular     <- lst_para$n_reuse_regular

    ## acceleration stage
    n_reuse_acc         <- lst_para$n_reuse_acc
    size_cohort_acc     <- lst_para$size_cohort_acc
    acc_max_dose        <- lst_para$acc_max_dose

    cur_dose            <- dose

    if (is.null(data)) {
        enroll_shift <- lst_para$date_bos
        sid_shift    <- 0
        n_treated    <- 0
        is_acc       <- TRUE
        pt_available <- data.frame()
    } else {
        enroll_shift <- max(data$date_dlt)
        sid_shift    <- max(data$sid)
        n_treated    <- nrow(data %>%
                             filter(dose == cur_dose))

        pt_level <- data %>%
            group_by(dose) %>%
            summarize(n = n())

        pt_treated <- data %>%
            filter(dose >= cur_dose) %>%
            select(sid) %>%
            unique()

        pt_available <- data %>%
            group_by(sid, date_enroll) %>%
            summarize(n    = n(),
                      ntox = sum(tox)) %>%
            filter(0 == ntox &
                   n < acc_max_dose) %>%
            filter(!(sid %in% pt_treated$sid)) %>%
            arrange(date_enroll)


        is_acc <- 0 == n_treated &
            max(pt_level$n) < size_cohort_regular
    }

    n_ava  <- nrow(pt_available)
    if (is_acc) {
        n_need  <- size_cohort_acc
        n_reuse <- n_reuse_acc
    } else {
        if (n_treated > 0 &
            n_treated < size_cohort_regular) {
            ## change from acc to regular
            n_need  <- size_cohort_regular - n_treated
        } else {
            n_need  <- size_cohort_regular
        }

        n_reuse <- n_reuse_regular
    }

    nenroll_exist <- min(n_need, n_reuse, n_ava)
    nenroll_new   <- n_need - nenroll_exist

    cur_new    <- NULL
    cur_exist  <- NULL

    cur_tox    <- lst_para$tox_rate[cur_dose]
    cur_res    <- lst_para$res_rate[cur_dose]
    rho        <- lst_para$rho

    if (nenroll_new > 0) {
        cur_prob <- stb_tl_p1_sce(
            c(cur_tox, cur_res, rho))[1, ]

        cur_new <- stb_tl_p1_simu_cohort(
            nenroll_new, cur_prob,
            par_enroll   = lst_para$par_enroll,
            enroll_shift = enroll_shift,
            dlt_days     = lst_para$dlt_days,
            ...)

        cur_new     <- cbind(cur_new$data_enroll,
                             cur_new$data_tox)
        cur_new$sid <- cur_new$sid + sid_shift
    }

    if (nenroll_exist > 0) {
        cur_prob <- stb_tl_p1_sce(
            c(cur_tox * lst_para$intra_fac, cur_res, rho))[1, ]

        cur_exist <- stb_tl_p1_simu_cohort(
            nenroll_exist, cur_prob,
            par_enroll   = lst_para$par_enroll,
            enroll_shift = 0,
            dlt_days     = lst_para$dlt_days,
            ...)

        sid_exist <- pt_available$sid[1:nenroll_exist]
        d2        <- cur_exist$data_tox %>%
            mutate(date_dlt = as.Date(enroll_shift + lst_para$dlt_days))

        d1_names  <- colnames(cur_exist$data_enroll)

        d1        <- data %>%
            filter(sid %in% sid_exist) %>%
            group_by(sid) %>%
            arrange(date_dlt, .by_group = TRUE) %>%
            slice_tail(n = 1) %>%
            select(all_of(d1_names))

        cur_exist <- cbind(d1, d2)
    }

    cur_data      <- rbind(cur_new, cur_exist)
    cur_data$dose <- dose

    cur_data
}


#' Generate Cohort
#'
#'
#'
desp1_generate_cohort_fix <- function(lst_para, dose, data, inx, ...) {
    cur_prob <- lst_para$theta[dose, ]
    cur_n    <- lst_para$size_cohort

    ## enrollment
    if (is.null(data)) {
        enroll_shift <- lst_para$date_bos
        sid_shift    <- 0
    } else {
        enroll_shift <- max(data$date_dlt)
        sid_shift    <- max(data$sid)
    }

    cur_data <- stb_tl_p1_simu_cohort(
        cur_n, cur_prob, inx,
        par_enroll   = lst_para$par_enroll,
        enroll_shift = enroll_shift,
        dlt_days     = lst_para$dlt_days,
        ...)

    cur_data$dose <- dose
    cur_data$sid  <- cur_data$sid + sid_shift

    ## return
    rbind(data, cur_data)
}

#' Summarize trial
#'
#' @export
#'
desp1_summarize <- function(data, n_dose, ...) {

    by_dose <- data %>%
        mutate(dose = factor(dose, levels = 1:n_dose)) %>%
        group_by(dose) %>%
        summarize(n_tot = n(),
                  n_tox = sum(tox),
                  n_res = sum(response)) %>%
        complete(dose,
                 fill = list(n_tot = 0, n_tox = 0, n_res = 0))

    by_study <- data.frame(
        n_tot    = length(unique(data$sid)),
        n_tox    = sum(data$tox),
        n_res    = sum(data$res),
        duration = as.Date(max(data$date_dlt)) -
            as.Date(min(data$date_enroll)),
        n_dose   = length(unique(data$dose)),
        n_cohort = length(unique(data$cohort)))

    list(by_dose  = by_dose,
         by_study = by_study)
}

#'  3+3 Design escalation
#'
#'
#' @export
#'
tp3_escalation <- function(lst_para, data, cur_dose) {

    n_dose    <- lst_para$n_dose
    cur_data  <- data %>% filter(dose == cur_dose)
    n_dlt     <- sum(cur_data$tox)
    n         <- nrow(cur_data)

    if (0 == n_dlt) {
        decision <- 1
    } else if (1 == n_dlt & 3 == n) {
        decision <- 0
    } else if (1 >= n_dlt) {
        decision <- 1
    } else {
        decision <- -1
    }

    next_dose <- min(cur_dose + decision, n_dose)
    next_dose <- max(next_dose, 1)

    tox_dose <- data %>%
        group_by(dose) %>%
        summarize(ntox = sum(tox)) %>%
        filter(ntox >= 2)

    if (next_dose %in% tox_dose$dose) {
        ## stop
        next_dose <- cur_dose
    }

    next_data <- data %>% filter(dose == next_dose)
    if (6 == nrow(next_data))
        next_dose <- -1

    next_dose
}

#' Accelerated design
#'
#' Generate cohorts
#'
acctit_generate_cohort <- function(lst_para, dose, data, inx, ...) {

    size_cohort  <- lst_para$size_cohort
    acc_max_dose <- lst_para$acc_max_dose
    cur_dose     <- dose

    pt_available <- data.frame()
    if (is.null(data)) {
        enroll_shift <- lst_para$date_bos
        sid_shift    <- 0
        n_need       <- 1
    } else {
        enroll_shift <- max(data$date_dlt)
        sid_shift    <- max(data$sid)

        pt_level <- data %>%
            group_by(dose) %>%
            summarize(n = n())

        pt_treated <- data %>%
            filter(dose >= cur_dose) %>%
            select(sid) %>%
            unique()

        pt_available <- data %>%
            group_by(sid, date_enroll) %>%
            summarize(n    = n(),
                      ntox = sum(tox)) %>%
            filter(0 == ntox &
                   n < acc_max_dose) %>%
            filter(!(sid %in% pt_treated$sid)) %>%
            arrange(date_enroll)

        n_treated <- nrow(data %>%
                          filter(dose == cur_dose))

        ## accelerated
        if (0 == n_treated &
            max(pt_level$n) < size_cohort) {
            n_need  <- 1
        } else {
            if (n_treated < size_cohort) {
                n_need <- size_cohort - n_treated
            } else {
                n_need <- size_cohort
            }
        }
    }

    n_ava   <- nrow(pt_available)
    n_new   <- max(n_need - n_ava, 0)
    n_exist <- n_need - n_new

    cur_new    <- NULL
    cur_exist  <- NULL
    cur_tox    <- lst_para$tox_rate[cur_dose]
    cur_res    <- lst_para$res_rate[cur_dose]
    rho        <- lst_para$rho

    if (n_new > 0) {
        cur_prob <- stb_tl_p1_sce(
            c(cur_tox, cur_res, rho))[1, ]

        cur_new <- stb_tl_p1_simu_cohort(
            n_new, cur_prob, inx,
            par_enroll   = lst_para$par_enroll,
            enroll_shift = enroll_shift,
            dlt_days     = lst_para$dlt_days,
            ...)
        cur_new$sid <- cur_new$sid + sid_shift
    }

    if (n_exist > 0) {
        cur_prob <- stb_tl_p1_sce(
            c(cur_tox * lst_para$intra_fac, cur_res, rho))[1, ]

        cur_exist <- stb_tl_p1_simu_cohort(
            n_exist, cur_prob, inx,
            par_enroll   = lst_para$par_enroll,
            enroll_shift = enroll_shift,
            dlt_days     = lst_para$dlt_days,
            ...)

        sid_exist <- pt_available$sid[1:n_exist]

        d1 <- data %>%
            filter(sid %in% sid_exist) %>%
            group_by(sid) %>%
            arrange(date_dlt, .by_group = TRUE) %>%
            slice_tail(n = 1) %>%
            mutate(date_dlt = date_dlt + lst_para$dlt_days) %>%
            select(1:6)

        d2 <- cur_exist %>%
            select(-(1:6))

        cur_exist <- cbind(d1, d2)
    }

    cur_data      <- rbind(cur_new, cur_exist)
    cur_data$dose <- dose

    ## return
    rbind(data, cur_data)
}


#'  3+3 Design escalation
#'
#'
#' @export
#'
acctit_escalation <- function(lst_para, data, cur_dose, ...) {

    size_cohort_acc <- lst_para$size_cohort_acc
    n_dose          <- lst_para$n_dose
    cur_data        <- data %>% filter(dose == cur_dose)
    n_dlt           <- sum(cur_data$tox)
    n               <- nrow(cur_data)

    if (size_cohort_acc == n & 0 == n_dlt) {
        next_dose <- min(cur_dose + 1, n_dose)
    } else if (size_cohort_acc == n) {
        next_dose <- cur_dose
    } else {
        next_dose <- tp3_escalation(lst_para, data, cur_dose, ...)
    }

    next_dose
}
