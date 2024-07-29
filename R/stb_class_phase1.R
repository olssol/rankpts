## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
##
##   DESCRIPTION:
##       DEFINE SIMULATION TOOLBOX CLASSES OF PHASE 1 DESIGNS
##
##
##   DESIGNS:
##     10: STB_DESIGN_DOSE_FIX
##
## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------

## -----------------------------------------------------------------------------
##                     10. Dose escaltion for Phase I
## -----------------------------------------------------------------------------
#'
#' @export
#'
setClass("STB_DESIGN_DOSE_P1",
         contains = "STB_DESIGN")

setMethod("stb_describe",
          "STB_DESIGN_DOSE_P1",
          function(x, ...) {
              desp1_describe(x, ...)
          })

setMethod("stb_set_default_para",
          "STB_DESIGN_DOSE_P1",
          function(x) {
              internal_desp1_dpara()
          })

setMethod("stb_para<-",
          "STB_DESIGN_DOSE_P1",
          function(x, value) {
              x     <- callNextMethod(x, value)
              theta <- stb_tl_p1_sce(x@design_para$tox_rate,
                                     x@design_para$res_rate,
                                     x@design_para$rho)
              x@design_para$theta <- theta
              x
          })

## -----------------------------------------------------------------------------
##                     10. Dose escaltion for FIX study
## -----------------------------------------------------------------------------
#'
#' @export
#'
setClass("STB_DESIGN_DOSE_FIX",
         contains = "STB_DESIGN")

setMethod("stb_describe",
          "STB_DESIGN_DOSE_FIX",
          function(x, ...) {
              callNextMethod()
              desfix_describe(x, ...)
          })

setMethod("stb_set_default_para",
          "STB_DESIGN_DOSE_FIX",
          function(x) {
              internal_desfix_dpara()
          })

setMethod("stb_plot_design",
          "STB_DESIGN_DOSE_FIX",
          function(x, ...) {
              desfix_plot_scenario(x@design_para, ...)
          })

setMethod("stb_generate_data",
          "STB_DESIGN_DOSE_FIX",
          function(x, ...) {
              desfix_gen_data(x@design_para, ...)
          })

setMethod("stb_analyze_data",
          "STB_DESIGN_DOSE_FIX",
          function(x, data_ana, ...) {
              rst <- desfix_single_trial(data_ana[[1]],
                                         lst_design = x@design_para,
                                         ...)

              list(rst)
          })

setMethod("stb_simu_gen_raw",
          "STB_DESIGN_DOSE_FIX",
          function(x, lst, ...) {
              n_reps  <- length(lst)
              rst     <- list()
              for (i in seq_len(length(lst))) {
                  rst[[i]]  <- lst[[i]][[1]] %>% mutate(rep = i)
              }

              list(n_reps = n_reps,
                   rst    = rbindlist(rst))
          })

setMethod("stb_simu_gen_summary",
          "STB_DESIGN_DOSE_FIX",
          function(x, lst, ...) {
              rst <- desfix_summary(
                  lst$rst, x@design_para, ...)
              list(rst)
          })
