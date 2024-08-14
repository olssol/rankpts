## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
##
##                                UI
##
## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
tip_txt <- function(label, message, placement = "bottom") {
    label |>
        span(
            tooltip(
                bs_icon("info-circle"),
                message,
                placement = placement))
}


panel_data <- list(
    card(
        card_body(fileInput("inDtaFile", "Choose Patients' Data File (.csv)",
                            accept = c(
                                "text/csv",
                                "text/comma-separated-values,text/plain",
                                ".csv")),
                  DT::dataTableOutput("dtData")),
        min_height  = 800,
        full_screen = TRUE
    )
)

panel_rules <- list(
    card(
        fileInput("inRuleFile", "Choose Ranking Rule File (.csv)",
                  accept = c(
                      "text/csv",
                      "text/comma-separated-values,text/plain",
                      ".csv")),
        DT::dataTableOutput("dtRule"),
        checkboxInput("inChkRule", "Check", FALSE),
        conditionalPanel(
            condition = "input.inChkRule == true",
            DT::dataTableOutput("dtRulePresent")
        )
    )
)

panel_ranked <- list(
    card(
        card_header("Ranked Patients"),
        DT::dataTableOutput("dtRanked")
    ),

    card(
        card_header("Selected Patients"),
        DT::dataTableOutput("dtRankedSel")
    )
)


panel_results <- list(
    layout_sidebar(
        sidebar = sidebar(
            width = 400,
            radioButtons("inSelVar",
                "Select Variable",
                "",
                selected = NULL)
        ),

        card(
            card_header("Observed Data"),
            plotOutput("plt_obs"),
            DT::dataTableOutput("dtRst")),

        height = 800
    )
)

tab_main <- function() {
    page_navbar(
        title   = "Rank Patients for Treatment Evaluation",
        theme   = bs_theme(
            bootswatch = "yeti",
            base_font  = font_google("Inter")
        ),

        nav_panel(title = "Load Data",   !!!panel_data),
        nav_panel(title = "Load Rules",  !!!panel_rules),
        nav_panel(title = "Ranked Data", !!!panel_ranked),
        nav_panel(title = "Results",     !!!panel_results)
    )
}


## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
##
##                                DATA
##
## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
get_pt_data <- reactive({
    req(input$inDtaFile)
    in_file <- input$inDtaFile
    read.csv(in_file$datapath)
})

observeEvent(input$inRuleFile, {
    in_file           <- input$inRuleFile
    userLog$data_rule <- read.csv(in_file$datapath)
})

proxy_dta_rule <- DT::dataTableProxy("dtRule")
observeEvent(input$dtRule_cell_edit, {
    info <- input$dtRule_cell_edit
    ## str(info)

    userLog$data_rule[info$row, info$col] <- info$value
    DT::replaceData(proxy_dta_rule,
                    userLog$data_rule,
                    resetPaging = FALSE)
})

get_rule_data <- reactive({
    req(userLog$data_rule)

    userLog$data_rule %>%
        filter(order > 0) %>%
        arrange(order)
})

get_cur_var <- reactive({
    cur_var <- input$inSelVar
    req(cur_var)

    if ("rank" == cur_var)
        return(cur_var)

    rule_dta <- get_rule_data()
    req(rule_dta)

    cur_var <- rule_dta %>%
        filter(variable == cur_var)

    cur_var
})

get_rank_test <- reactive({
    pt_dta <- get_pt_data()
    req(pt_dta)


    rule_dta <- get_rule_data()
    req(rule_dta)

    rpt_rank_test(pt_dta, rule_dta)
})

get_rank_test_data <- reactive({
    rst <- get_rank_test()
    req(rst)
    rst$data
})

get_rank_test_data_selected <- reactive({
    ranked_data <- get_rank_test_data()
    req(ranked_data)

    selected <- input$dtRanked_rows_selected
    req(selected)

    rule_dta <- get_rule_data()
    req(rule_dta)

    rpt_get_rank_data(ranked_data, rule_dta, selected)
})


get_rank_test_tbl <- reactive({
    rst <- get_rank_test()
    req(rst)
    rst$test
})

get_rank_test_plot <- reactive({
    rst <- get_rank_test_data()
    req(rst)

    rpt_plot_hist(rst, vname = "rank", vlabel = "Rank")
})


get_test_rst <- reactive({
    pt_dta <- get_pt_data()
    req(pt_dta)

    cur_var <- get_cur_var()
    req(cur_var)

    if (!is.data.frame(cur_var)) {
        if ("rank" == cur_var) {
            return(get_rank_test_tbl())
        }
    }

    if ("continuous" == cur_var[1, "type"]) {
        ftest <- rpt_ana_ttest
    } else if ("binary" == cur_var[1, "type"]){
        ftest <- rpt_ana_fisher
    } else {
        return(NULL)
    }

    pt_dta$y <- pt_dta[[cur_var[1, "variable"]]]
    rst      <- ftest(pt_dta)
    row.names(rst) <- NULL
    rst
})

plot_rst <- reactive({
    pt_dta <- get_pt_data()
    req(pt_dta)

    cur_var <- get_cur_var()
    req(cur_var)


    if (!is.data.frame(cur_var)) {
        return(get_rank_test_plot())
    }

    rpt_plot_hist(pt_dta, vname = cur_var[1, "variable"],
                  vlabel = cur_var[1, "label"])

})
