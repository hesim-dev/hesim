# The UI.R script contains the ui, this is the layout of the graphical user interface (GUI) of your app
# This can be read linearly, going through content positioning tab by tab

# Load packages ---------------
# Packages will need installing if they have not been installed before (use install.packages())
# Packages can be defined in either the server.R or ui.R scripts, but ui.R is read first
library(shiny)             # Shiny package to produce the app
library(shinydashboard)    # Functions to give the app a tabbed layout
library(hesim)             # Required for the model functionality (see the 'hesim example model.R script')
library(data.table)        # Required for the model functionality (see the 'hesim example model.R script')
library(flexsurv)          # Required for the model functionality (see the 'hesim example model.R script')
library(ggplot2)           # Required for creating graphs
library(magrittr)          # Required for the model functionality (see the 'hesim example model.R script')
library(shinyWidgets)      # Used for creating lovely input options (in this case the % on the discount)
                           # see http://shinyapps.dreamrs.fr/shinyWidgets/ to see examples
library(diagram)           # Assists with creating the diagram from heemod
library(heemod)            # Can produce a really simple model diagram - also has other useful functions for partitioned survival modelling
library(survminer)         # useful for easily presenting Kaplan–Meier plots.
library(shinycssloaders)   # Creates the loading animation
library(DT)                # Creates datatables
library(shinyjs)           # Adds functionality using a library of javascript functions. E.g. click() can be used to click
                           # buttons in the back end without user interaction. This is user in server.R to start the model calculations running

options(encoding = "UTF-8")

ui <-
  #This app uses the shinydashboard package to create a layout - see https://rstudio.github.io/shinydashboard/ for details
  dashboardPage(
    # ~ Header -------------------
    dashboardHeader(
      title = "",
      titleWidth = 0 ,
      tags$li(
        class = "dropdown",
        style = " font-size: 16px; font-weight: bold;",
        tags$a(("IndivCtstm demo shiny app"))
      )
    ),
    
    # ~ Sidebar -------------------
    # This is not used in this model version, but can be added and laid out via the dashboardSidebar() function
    dashboardSidebar(disable = TRUE),
    
    # ~ Body -------------------
    # Arranging the content and outputs within the headings and sub-headings
    dashboardBody(
      shinyjs::useShinyjs(), #this is needed to use the click() functions, and other functions used for enabling/disabling inputs
      # ~~ Switches -------------------
      fluidRow(
        box(width = 3, 
            title = "Switches", 
            status = "primary", 
            solidHeader = TRUE,
            collapsible = TRUE,
            sliderInput(inputId = "Input_nSamples", 
                        label = "Number of probabilistic samples:",
                        min = 100, 
                        max = 1000, 
                        value = 100),
            sliderInput(inputId = "Input_timehorizon", 
                        label = "Enter time horizon (years):",
                        min = 10, 
                        max = 30, 
                        value = 30),
            numericInputIcon(
              inputId = "Input_discount_QALY",
              label = "Discount for QALYs:",
              min = 0,
              max = 100,
              value = 1,
              icon = list(NULL, icon("percent"))
            ),
            numericInputIcon(
              inputId = "Input_discount_Costs",
              label = "Discount for Costs:",
              min = 0,
              max = 100,
              value = 3.5,
              icon = list(NULL, icon("percent"))
            ),
            column(12,hr()),
            column(12,
                   actionButton("Run_model",
                                "Re-run the model", icon("sync"),
                                style = "color: #ffffff; background-color: #222D32; border-color: #1C75BB"), align = 'center'),
            column(12,br()),
            column(12,
                   downloadButton("Create_htmlreport",
                                  "Create model html report", 
                                  style = "color: #ffffff; background-color: #222D32; border-color: #1C75BB"), align = 'center'),
            column(12,br()),
            column(12,
                   downloadButton("Create_pdfreport",
                                  "Create model pdf report", 
                                  style = "color: #ffffff; background-color: #222D32; border-color: #1C75BB"), align = 'center')
        ),
        # ~~ Inputs -------------------
        tabBox(
          title = "Inputs",
          id = "Input_tabBox",
          selected = "States and transitions",
          width = 9,side = "right",
          tabPanel(title = "Utilities",
                   # ANSWER CODE CHANGES ARE HERE ----------------------------
                   fluidRow(
                     column(12,"Text can be written here to explain how the inputs are used if the coder chooses",br(),br()),
                     column(6, tags$u("Utility inputs"),br(),tableOutput("Utility_out"))
                   )),
          tabPanel(title = "Costs",
                   fluidRow(
                     column(8, tags$u("Drug costs per strategy, state and time"),br(),tableOutput("Drugcost_out")),
                     column(4, tags$u("Medical costs per state"),br(),tableOutput("Medcost_out"))
                   )),
          tabPanel(title = "Patients",
                   fluidRow(
                     column(12,
                            htmlOutput("Patient_number"),
                            plotOutput("Patient_hist"))
                   )),
          tabPanel(title = "Trial data",
                   fluidRow(
                     column(12,
                            selectInput("trial_trans_input",
                                        "Select transition to view data",
                                        choices = c(1,2,3)),
                            withSpinner(# this is from the shinycssloaders package
                              plotOutput("Trial_data_plot", height = "500px")))
                   )),
          tabPanel(title = "Strategies",
                   fluidRow(
                     column(6, tags$u("Model strategies"),br(),tableOutput("strategies_out"))
                   )),
          tabPanel(title = "Model diagram",
                   fluidRow(
                     align = "center",
                     plotOutput(
                       outputId = "diagram",
                       width = "620px",
                       height = "450px"
                     )
                   )), 
          tabPanel(title = "States and transitions",
                   fluidRow(
                     column(6, tags$u("State names"),br(),tableOutput("state_out")),
                   column(6, tags$u("Transition matrix"),br(),tableOutput("tmat_out")),
                   column(6, tags$u("Transition IDs"),br(),tableOutput("transitions_out"))
                   ))
        ),
        # ~~ Results -------------------
        box(
          title = "Results",
          width = 12,
          status = "primary",
          solidHeader = TRUE,
          collapsible = TRUE,
          collapsed = FALSE,
          fluidRow(column(
            12,
            htmlOutput("Results_text"), br(),br(),
            HTML("<u><b>Survival outcomes </b></u>"),
            br(),
            withSpinner(plotOutput("Results_graph")),
            br(),
            br(),
            HTML("<u><b>Summary table </b></u>"),
            br(),
            withSpinner(dataTableOutput("Results_DT"))
          ))
        )
      )
    )
  )