#' TTest
#'
#' @export
#'
rpt_ana_ttest <- function(dta, varm = "ARM", vtrt = "Treatment") {
    d_trt <- dta %>% filter(!!sym(varm) == vtrt)
    d_ctl <- dta %>% filter(!!sym(varm) != vtrt)

    rst_test <- t.test(d_trt$y, d_ctl$y)

    data.frame(
        Test      = "T Test",
        Mean_Trt  = rst_test$estimate[1],
        Mean_Ctl  = rst_test$estimate[2],
        Conf_Low  = rst_test$conf.int[1],
        Conf_High = rst_test$conf.int[2],
        p_Value   = rst_test$p.value)
}

#' TTest
#'
#' @export
#'
rpt_ana_fisher <- function(dta, varm = "ARM", vtrt = "Treatment") {
    d_trt <- dta %>% filter(!!sym(varm) == vtrt)
    d_ctl <- dta %>% filter(!!sym(varm) != vtrt)

    cont_table <- table(dta$y, dta[[varm]])
    rst_test   <- fisher.test(cont_table)

    data.frame(
        Test      = "Fisher Exact",
        Mean_Trt  = mean(d_trt$y, na.rm = TRUE),
        Mean_Ctl  = mean(d_ctl$y, na.rm = TRUE),
        Conf_Low  = rst_test$conf.int[1],
        Conf_High = rst_test$conf.int[2],
        p_Value   = rst_test$p.value)
}


#' TTest
#'
#' @export
#'
rpt_plot_hist <- function(dta, vname, vlabel, varm = "ARM", ...) {
    ggplot(dta, aes_string(x = vname, fill = varm)) +
        geom_histogram(
            aes(y    = after_stat(density)),
            position = "dodge",
            color    = "black",
            ...) +
        geom_density(size = 1, alpha = 0.4) +
        labs(x = vlabel,
             y = "Frequency") +
        theme_bw()
}
