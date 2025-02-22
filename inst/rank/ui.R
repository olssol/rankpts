library(shiny)
library(bslib)
library(ggplot2)
library(bsicons)
library(dplyr)
library(rankpts)

shinyUI(
    div(
        tags$head(tags$title("Rank Patients"),
                  tags$link(rel = "stylesheet", type = "text/css",
                            href = "styles.css")
                  ),

        ##wait for starting
        conditionalPanel(condition = "$('html').hasClass('shiny-busy')",
                         tags$div("", id = "loadmessage")),

        ##main page
        uiOutput("mainpage"),

        ##foot
        withTags({
            div(class = "cfooter",
                "A",
                a("Statistical Innovation", href = "http://www.regeneron.com/"),
                "Project")
        })
    )
)
