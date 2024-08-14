shinyServer(function(input, output, session) {

    source("design_ui.R", local = TRUE);

    userLog           <- reactiveValues()
    userLog$data_rule <- NULL

    ##---------main-------------------------
    ##--------------------------------------
    output$mainpage <- renderUI({
        tab_main()
    })

    ##---------patient data-----------------
    ##--------------------------------------
    output$dtData <- DT::renderDataTable({
                             pt_data <- get_pt_data()
                             req(pt_data)
                             pt_data
                         },
                         options = list(pageLength = 100)
                         )

    ##---------rank rules data--------------
    ##--------------------------------------
    output$dtRule <- DT::renderDataTable({
                             rule_data <- userLog$data_rule
                             req(rule_data)
                             rule_data
                         },
                         editable = list(target  = 'cell',
                                         disable = list(columns = c(1, 2))),
                         server   = TRUE,
                         options  = list(
                             ordering = FALSE,
                             paging   = FALSE,
                             dom      = 't',
                             info     = FALSE)
                         )

    output$dtRulePresent <- DT::renderDataTable({
                             rule_data <- get_rule_data()
                             req(rule_data)
                             rule_data
                         },
                         options  = list(
                             paging   = FALSE,
                             dom      = 't',
                             info     = FALSE)
                         )

    ##---------ranked data-----------------
    ##--------------------------------------
    output$dtRanked <- DT::renderDataTable({
                               ranked <- get_rank_test_data()
                               req(ranked)
                               ranked
                           },
                           selection = 'multiple',
                           options = list(pageLength = 100)
                           )

    output$dtRankedSel <- DT::renderDataTable({
                                    data <- get_rank_test_data_selected()
                                    req(data)
                                    data
                                },
                                options  = list(
                                    paging   = FALSE,
                                    dom      = 't',
                                    info     = FALSE)
                                )

    ##---------result-----------------------
    ##--------------------------------------
    observeEvent(get_rule_data(), {
        cur_sel  <- input$inSelVar
        rule_dta <- get_rule_data()
        req(rule_dta)

        vars    <- rule_dta$variable
        labels  <- rule_dta$label

        vars    <- c(vars, "rank")
        labels  <- c(labels, "RANK BASED ON ALL VARIABLES")

        names(vars) <- labels
        updateRadioButtons(session,
                           "inSelVar",
                           choices  = vars,
                           selected = cur_sel)
    })

    output$plt_obs <- renderPlot({
        plot_rst()
    })

    output$dtRst <- DT::renderDataTable({
                            rst <- get_test_rst()
                            req(rst)
                            rst
                        },
                        options  = list(
                            paging   = FALSE,
                            dom      = 't',
                            info     = FALSE)
                        )

})
