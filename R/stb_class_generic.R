## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
##
##   DESCRIPTION:
##       DEFINE SIMULATION TOOLBOX GENERIC FUNCTIONS
##
##
## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------


## -----------------------------------------------------------------------------
##                  class stb_design
## -----------------------------------------------------------------------------

setGeneric("stb_describe",
           function(x, ...) standardGeneric("stb_describe"))

setGeneric("stb_get_para",
           function(x) standardGeneric("stb_get_para"))

setGeneric("stb_para<-",
           function(x, value) standardGeneric("stb_para<-"))

setGeneric("stb_set_default_para",
           function(x, ...) standardGeneric("stb_set_default_para"))

setGeneric("stb_plot_design",
           function(x, ...) standardGeneric("stb_plot_design"))

setGeneric("stb_generate_data",
           function(x, ...) standardGeneric("stb_generate_data"))

setGeneric("stb_generate_cohort",
           function(x, dose, data, ...) standardGeneric("stb_generate_cohort"))

setGeneric("stb_escalation",
           function(x, dose, data, ...) standardGeneric("stb_escalation"))

setGeneric("stb_recommend",
           function(x, data, ...) standardGeneric("stb_recommend"))

setGeneric("stb_analyze_data",
           function(x, data, ...) standardGeneric("stb_analyze_data"))

setGeneric("stb_plot_data",
           function(x, data, ...) standardGeneric("stb_plot_data"))

setGeneric("stb_create_trial",
           function(x, ...) standardGeneric("stb_create_trial"))

setGeneric("stb_create_analysis_set",
           function(x, data, ...) standardGeneric("stb_create_analysis_set"))

setGeneric("stb_create_simustudy",
           function(x,
                    n_rep  = 5,
                    n_core = 5,
                    seed   = NULL, ...) standardGeneric("stb_create_simustudy"))

setGeneric("stb_simu_gen_raw",
           function(x, lst, ...) standardGeneric("stb_simu_gen_raw"))

setGeneric("stb_simu_gen_summary",
           function(x, lst, ...) standardGeneric("stb_simu_gen_summary"))

setGeneric("stb_simu_gen_key",
           function(x, lst, ...) standardGeneric("stb_simu_gen_key"))


## -----------------------------------------------------------------------------
##                        class stb_trial
## -----------------------------------------------------------------------------

setGeneric("stb_get_trial_data",
           function(x) standardGeneric("stb_get_trial_data"))

setGeneric("stb_get_trial_result",
           function(x) standardGeneric("stb_get_trial_result"))

setGeneric("stb_get_trial_seed",
           function(x) standardGeneric("stb_get_trial_seed"))

setGeneric("stb_get_trial_design",
           function(x) standardGeneric("stb_get_trial_design"))

setGeneric("stb_trial_plot",
           function(x, ...) standardGeneric("stb_trial_plot"))


## -----------------------------------------------------------------------------
##                        class stb_simustudy
## -----------------------------------------------------------------------------

setGeneric("stb_get_simu_design",
           function(x) standardGeneric("stb_get_simu_design"))

setGeneric("stb_get_simu_raw",
           function(x) standardGeneric("stb_get_simu_raw"))

setGeneric("stb_get_simu_summary",
           function(x) standardGeneric("stb_get_simu_summary"))

setGeneric("stb_get_simu_key",
           function(x) standardGeneric("stb_get_simu_key"))

setGeneric("stb_get_simu_nrep",
           function(x) standardGeneric("stb_get_simu_nrep"))

setGeneric("stb_get_simu_seed",
           function(x) standardGeneric("stb_get_simu_seed"))



## -----------------------------------------------------------------------------
##                        overall class stb_design
## -----------------------------------------------------------------------------

#'
#' @export
#'
setClass("STB_DESIGN",
         slots     = list(design_para  = "list",
                          design_valid = "numeric"),
         prototype = prototype(design_para  = list(),
                               design_valid = 1))

setMethod("initialize",
          "STB_DESIGN",
          function(.Object, ...) {
              .Object@design_para <- stb_set_default_para(.Object)
              .Object
          })

#'
#' @export
#'
setMethod("stb_set_default_para",
          "STB_DESIGN",
          function(x) list())

#'
#' @export
#'
setMethod("stb_get_para",
          "STB_DESIGN",
          function(x) x@design_para)

#'
#' @export
#'
setMethod("stb_para<-",
          "STB_DESIGN",
          function(x, value) {
              x@design_para <- tl_merge_lists(value,
                                              x@design_para)
              x
          })

#'
#' @export
#'
setMethod("stb_describe",
          "STB_DESIGN",
          function(x, ...) NULL)

#'
#' @export
#'
setMethod("stb_generate_data",
          "STB_DESIGN",
          function(x, ...) NULL)

#'
#' @export
#'
setMethod("stb_generate_cohort",
          "STB_DESIGN",
          function(x, dose, data, ...) NULL)

#'
#' @export
#'
setMethod("stb_escalation",
          "STB_DESIGN",
          function(x, dose, data, ...) NULL)

#'
#' @export
#'
setMethod("stb_recommend",
          "STB_DESIGN",
          function(x, data, ...) NULL)

#'
#' @export
#'
setMethod("stb_analyze_data",
          "STB_DESIGN",
          function(x, data, ...) list())

#'
#' @export
#'
setMethod("stb_create_trial",
          "STB_DESIGN",
          function(x, seed = NULL, ...) {

    if (!is.null(seed))
        old_seed <- set.seed(seed)

    next_dose <- 1
    data      <- NULL
    inx       <- 1
    while (next_dose > -1) {
        dose      <- next_dose
        data      <- stb_generate_cohort(x, dose, data, inx, ...)
        next_dose <- stb_escalation(x, dose, data, ...)
        inx       <- inx + 1
    }

    result                    <- stb_analyze_data(x, data, ...)
    result$by_study$recommend <- stb_recommend(x, data)

    if (!is.null(seed))
        set.seed(old_seed)

    new("STB_TRIAL",
        design    = x,
        data      = data,
        result    = result,
        seed      = seed)
})

#'
#' @export
#'
setMethod("stb_simu_gen_raw",
          "STB_DESIGN",
          function(x, lst, ...) {

    by_dose  <- list()
    by_study <- list()

    for (i in seq_len(length(lst))) {
        cur_rst       <- lst[[i]]
        by_dose[[i]]  <- cur_rst$by_dose %>% mutate(rep = i)
        by_study[[i]] <- cur_rst$by_study %>% mutate(rep = i)
    }

    list(by_dose  = rbindlist(by_dose),
         by_study = rbindlist(by_study))
})

#'
#' @export
#'
setMethod("stb_simu_gen_summary",
          "STB_DESIGN",
          function(x, lst, ...) {

    ## by dose
    by_dose <- lst$by_dose %>%
        group_by(dose) %>%
        summarize(n_tot = mean(n_tot),
                  n_tox = mean(n_tox),
                  n_res = mean(n_res))

    ## by study
    by_study <- lst$by_study %>%
        summarize(n_tot    = mean(n_tot),
                  n_tox    = mean(n_tox),
                  n_res    = mean(n_res),
                  duration = mean(duration),
                  n_dose   = mean(n_dose),
                  n_cohort = mean(n_cohort))

    ## recommend
    dose_recommend <- lst$by_study %>%
        group_by(recommend) %>%
        summarize(freq = n()) %>%
        mutate(freq = freq / sum(freq))

    list(by_dose        = by_dose,
         by_study       = by_study,
         dose_recommend = dose_recommend)
})

#'
#' @export
#'
setMethod("stb_simu_gen_key",
          "STB_DESIGN",
          function(x, lst, ...) lst)


#'
#' @export
#'
setMethod("stb_create_simustudy",
          "STB_DESIGN",
          function(x,
                   n_rep    = 5,
                   n_core   = 5,
                   seed     = NULL,
                   ...,
                   save_raw = FALSE) {

              if (0 == x@design_valid) {
                  cat("Invalid design. \n")
                  return(NULL)
              }

              ## seed
              if (!is.null(seed))
                  old_seed <- set.seed(seed)

              ## all random seeds
              all_seeds <- ceiling(abs(rnorm(n_rep) * 100000))

              ## replications
              rst <-
                  parallel::mclapply(
                                seq_len(n_rep),
                                function(k) {
                                    if (0 == k %% 5)
                                        print(k)

                                    cur_trial <-
                                        stb_create_trial(
                                            x,
                                            seed = all_seeds[k],
                                            ...)

                                    cur_trial@result
                                },
                                mc.cores = n_core)

              ## summarize
              rst_raw     <- stb_simu_gen_raw(x, rst)
              rst_summary <- stb_simu_gen_summary(x, rst_raw, ...)

              ## seed
              if (!is.null(seed))
                  set.seed(old_seed)

              if (!save_raw)
                  rst_raw <- list()

              new("STB_SIMU_STUDY",
                  design      = x,
                  n_rep       = n_rep,
                  n_core      = n_core,
                  rst_raw     = rst_raw,
                  rst_summary = rst_summary,
                  seed        = seed)
          })


## -----------------------------------------------------------------------------
##                        overall class stb_trial
## -----------------------------------------------------------------------------

#'
#'
#' @export
#'
setClass("STB_TRIAL",
         slots = list(design      = "STB_DESIGN",
                      data        = "data.frame",
                      result      = "list",
                      seed        = "numeric"),
         prototype = prototype(seed     = NULL,
                               result   = NULL))

#'
#' @export
#'
setMethod("stb_get_trial_data",   "STB_TRIAL", function(x) x@data)

#'
#' @export
#'
setMethod("stb_get_trial_result", "STB_TRIAL", function(x) x@result)

#'
#' @export
#'
setMethod("stb_get_trial_design", "STB_TRIAL", function(x) x@design)

#'
#' @export
#'
setMethod("stb_get_trial_seed",   "STB_TRIAL", function(x) x@seed)


#'
#' @export
#'
setMethod("stb_trial_plot",
          "STB_TRIAL",
          function(x, ...) {
              stb_plot_data(x@design, x@data, ...)
          })


## -----------------------------------------------------------------------------
##                        overall class stb_simustudy
## -----------------------------------------------------------------------------

#'
#' @export
#'
setClass("STB_SIMU_STUDY",
         slots = list(design      = "STB_DESIGN",
                      n_rep       = "numeric",
                      n_core      = "numeric",
                      rst_raw     = "list",
                      rst_summary = "list",
                      rst_key     = "list",
                      seed        = "numeric"))


#'
#' @export
#'
setMethod("stb_get_simu_design",   "STB_SIMU_STUDY", function(x) x@design)

#'
#' @export
#'
setMethod("stb_get_simu_raw",      "STB_SIMU_STUDY", function(x) x@rst_raw)

#'
#' @export
#'
setMethod("stb_get_simu_summary",  "STB_SIMU_STUDY", function(x) x@rst_summary)

#'
#' @export
#'
setMethod("stb_get_simu_key",      "STB_SIMU_STUDY", function(x) x@rst_key)

#'
#' @export
#'
setMethod("stb_get_simu_nrep",     "STB_SIMU_STUDY", function(x) x@n_rep)

#'
#' @export
#'
setMethod("stb_get_simu_seed",     "STB_SIMU_STUDY", function(x) x@seed)



## -----------------------------------------------------------------------------
##       HELPER FUNCTIONS
## -----------------------------------------------------------------------------

#' Create a Study Design
#'
#'
#'
#' @export
#'
stb_create_design <- function(type = c("dose_fix",
                                       "dose_3p3",
                                       "dose_acctit",
                                       "dose_acc2")) {

    type <- match.arg(type)
    rst  <- switch(type,
                   dose_fix    = new("STB_DESIGN_DOSE_FIX"),
                   dose_3p3    = new("STB_DESIGN_P1_3P3"),
                   dose_acctit = new("STB_DESIGN_P1_ACCTIT"),
                   dose_acc2   = new("STB_DESIGN_P1_ACC2"),
                   new("STB_DESIGN"))

    rst
}
