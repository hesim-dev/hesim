#The server.R script is the model functionality. It is laid out in the same way as the 'hesim example model.R' script
# Unlike the 'hesim example model.R' script, the server takes in inputs from the 'inputs' defined in the ui.R, and used these
# in the calculations. It then produced the 'outputs' as defined in the ui.R, these are rendered here to appear in the 
# shiny model graphical user interface

library("shiny")
options(encoding = "UTF-8")


server <- function(input, output, session) { 
  # Informing inputs --------------
  # ~ states and transitions ---------
  # ~~ Define matrix ---------------
  tmat <- rbind(
    c(NA, 1, 2),
    c(NA, NA, 3),
    c(NA, NA, NA)
  )
  colnames(tmat) <- rownames(tmat) <- c("Stable", "Progression", "Death")
  #print(tmat)
  
  # ~~ Define transitions ---------------
  transitions <- create_trans_dt(tmat)
  
  # ~~ Outline states and IDs in separate table for easy referencing --------------
  # Death is automatically added by get_labels() (below) in the code below in the default settings,
  # but 'death_label = NULL' argument in get_labels() this will override this. Current setup is to maintain simplicity
  states <- data.table(
    state_id = 1:2,
    state_name = c("Stable", "Progression")
  )
  
  #I have added an extra table here to present in output$state_out. Remember, the front-end user cannot see the back-end
  # comments, so as much information needs to be made available front end as possible
  states_wDeath <- data.table(
    state_id = 1:3,
    state_name = c("Stable", "Progression", "Death")
  )
  
  # ~~ Output to shiny -----------
  output$state_out <- renderTable(
    states_wDeath,
    striped = TRUE,
    digits = 0,
    bordered = TRUE,
    colnames = TRUE,
    sanitize.text.function = function(x) x) 
  
  output$tmat_out <- renderTable(
    tmat,
    striped = TRUE,
    digits = 0,
    bordered = TRUE,
    colnames = TRUE,
    rownames = TRUE,
    sanitize.text.function = function(x) x)
  
  output$transitions_out <- renderTable(
    transitions,
    striped = TRUE,
    digits = 0,
    bordered = TRUE,
    colnames = TRUE,
    sanitize.text.function = function(x) x) 
  
  Model_Diagram <- define_transition( #this function is part of the heemod package
    state_names = c("Stable", "Progression", "Death"),
    Stable,transition_id_1, transition_id_2,
    ,Progressed, transition_id_3,
    , , Death
  )
  
  output$diagram <- renderPlot({
    plot(Model_Diagram)
  })  
  
  
  # ~ Strategies ----------------------
  # ~~ Outline strategy and IDs ----------------
  strategies <- data.table(
    strategy_id = 1:3,
    strategy_name = c("SOC", "New 1", "New 2")
  )
  
  # ~~ Output to shiny -------------
  
  output$strategies_out <- renderTable(
    strategies,
    striped = TRUE,
    digits = 0,
    bordered = TRUE,
    colnames = TRUE,
    sanitize.text.function = function(x) x) 
  
  # ~ Patients -------------
  # ~~ Create patient sample to model -------------
  n_patients <- 1000
  patients <- data.table(
    patient_id = 1:n_patients,
    age = rnorm(n_patients, mean = 45, sd = 7),
    female = rbinom(n_patients, size = 1, prob = .51)
  )
  # If groups are wanted, these can be defined in the 'grp_id' and 'grp_name' columns. Otherwise can be commented and left blank.
  # patients[, grp_id := ifelse(female == 1, 1, 2)]
  # patients[, grp_name := ifelse(female == 1, "Female", "Male")]
  
  
  # ~~ Output to shiny -------------
  
  output$Patient_number <- renderText(
    HTML("The number of patients simulated for this model is: <b> ",n_patients,"</b>")
  )
  
  output$Patient_hist <- renderPlot(
    ggplot(patients, aes(x = age, fill = as.factor(female))) + 
      geom_histogram(binwidth = 1, colour = "#959595") + 
      theme_bw() + 
      scale_fill_manual("Gender:", values = c("#0D8E1E","#9552BB"), labels = c("Male","Female"))
    )
  
  # ~ Organising basic model settings ------------  
  # ~~ Create hesim data object -----------
  hesim_dat <- hesim_data(
    strategies = strategies,
    patients = patients,
    states = states,
    transitions = transitions
  )
  
  #print(hesim_dat)
  
  # ~~ Setting up labels for state and strategy IDs ---------------
  labs_indiv <- get_labels(hesim_dat)
  #print(labs_indiv)
  
  # ~ 'Trial' data ----------
  # hesim package includes the 'onc3' data.table. This separates the three transitions by 'transition_id', where the IDs match the 'transitions' data
  # These individual transitions can be filtered for and have parametric models fitted
  # Data example showing patients 1 and 2:
  # onc3[patient_id %in% c(1, 2)]
  
  # ~~ Fit the survival data ------------------
  n_trans <- max(tmat, na.rm = TRUE)
  wei_fits <- vector(length = n_trans, mode = "list")
  f <- as.formula(Surv(time, status) ~ factor(strategy_name) + female + age)
  
  for (i in 1:length(wei_fits)){
    if (i == 3) {f <- update(f, .~.-factor(strategy_name))} 
    wei_fits[[i]] <- flexsurvreg(f, data = onc3,
                                 subset = (transition_id == i),
                                 dist = "weibull")
  }
  
  wei_fits <- flexsurvreg_list(wei_fits)
  
  # ~~ Output to shiny -------------
  output$Trial_data_plot <- renderPlot({
    req(input$trial_trans_input)
    
    transition_id_view <- input$trial_trans_input
    TransitionData <-  survfit(as.formula(Surv(time, status) ~ strategy_name), data = onc3[which(transition_id == transition_id_view), ])
    ggsurvplot(
      fit      = TransitionData,
      data     = onc3,
      # break.y.by = 0.1,
      # break.x.by = 0.5,
      xlab = 'Time (Years)',
      #xlim = c(0,5),
      ylab = 'Survival',
      #palette=c("red","blue","green"),
      risk.table = TRUE,
      #risk.table.y.text.col = TRUE,
      #risk.table.height = 0.3,
      #risk.table.title = 'Number at risk',
      #conf.int = T,
      #linetype = c(1,2),
      legend = "top"
    )
  })
  
  # ~ Costs ---------------------
  # ~~ Create time-dependent drug costs per strategy --------------
  # The time units are in years.
  drugcost_dt <- matrix(c(
    1, 1, 1, 0.00, 0.25,  2000,
    1, 1, 2, 0.25, Inf, 2000,
    1, 2, 1, 0.00, 0.25,  1500,
    1, 2, 2, 0.25, Inf , 1200,
    2, 1, 1, 0.00, 0.25,  12000,
    2, 1, 2, 0.25, Inf , 12000,
    2, 2, 1, 0.00, 0.25,  1500,
    2, 2, 2, 0.25, Inf , 1200,
    3, 1, 1, 0.00, 0.25,  15000,
    3, 1, 2, 0.25, Inf , 15000,
    3, 2, 1, 0.00, 0.25,  1500,
    3, 2, 2, 0.25, Inf , 1200
  ),byrow = TRUE, ncol = 6, dimnames = list(NULL, c("strategy_id", "state_id", "time_id", "time_start", "time_stop","est")))
  drugcost_dt <- data.table(drugcost_dt)
  #print(drugcost_dt)
  
  drugcost_tbl <- stateval_tbl(
    drugcost_dt,
    dist = "fixed")
  #print(drugcost_tbl)
  
  # ~~ Medical costs ---------------
  medcost_tbl <- stateval_tbl(
    data.table(state_id = states$state_id,
               mean = c(2000, 9500),
               se = c(2000, 9500)
    ),
    dist = "gamma")
  #print(medcost_tbl)
  
  # ~~ Output to shiny -------------
  
  output$Drugcost_out <- renderTable(
    drugcost_tbl,
    striped = TRUE,
    bordered = TRUE,
    colnames = TRUE,
    sanitize.text.function = function(x) x)
  
  output$Medcost_out <- renderTable(
    medcost_tbl,
    striped = TRUE,
    bordered = TRUE,
    colnames = TRUE,
    sanitize.text.function = function(x) x)
  
  # ~ Utilities ----------------
  utility_tbl <- stateval_tbl(
    data.table(state_id = states$state_id,
               mean = c(.8, .6),
               se = c(0.02, .05)
    ),
    dist = "beta")
  
  # ~~ Output to shiny -------------
  # ANSWER CODE CHANGES ARE HERE ----------------------------
  output$Utility_out <- renderTable(
    utility_tbl,
    striped = TRUE,
    bordered = TRUE,
    colnames = TRUE,
    sanitize.text.function = function(x) x)
  #print(utility_tbl)
  
  
  # Setting up the model --------------
  # ~ Taking in inputs from shiny interface ---------------------------
  # ~~ Number of parameter samples is needed to use for the PSA
  
  # All of these elements are reactive, so they change whenever they are interacted with. You can view the 
  # value changing and printing to the R console with an observe() (see below)
  n_samples <- reactive({
    req(input$Input_nSamples)
    input$Input_nSamples
  })
  
  # You can have an observe to print input$Input_nSamples to the console as feedback
  #observe({print(paste0("The current number of samples is: ", input$Input_nSamples))})
  
  # ~~ Years in time horizon
  n_years <- reactive({
    req(input$Input_timehorizon)
    input$Input_timehorizon
  })
  
  # ~~ Discount QALY
  disc_QALY <- reactive({
    req(input$Input_discount_QALY)
    input$Input_discount_QALY/100
  })
  
  # ~~ Discount costs
  disc_Cost <- reactive({
    req(input$Input_discount_Costs)
    input$Input_discount_Costs/100
  })
  
  
  # ~ Expanding the data input dataframe to set up for running all patients with all strategies
  transmod_data <- expand(hesim_dat,
                          by = c("strategies", "patients"))
  head(transmod_data)
  
  # ~ Wrapping inputs in hesim functions for use in model -------------
  # ~~ Efficacy -----------------
  # these objects contain reactive data (e.g. n_samples()), so therefore need to be reactive themselves
  transmod <- reactive(create_IndivCtstmTrans(wei_fits, transmod_data,
                                     trans_mat = tmat, n = n_samples(),
                                     uncertainty = "normal",
                                     clock = "reset",
                                     start_age = patients$age))
  
  # ~~ Utilities -----------------
  utilitymod <- reactive(create_StateVals(utility_tbl, n = n_samples(),
                                 hesim_data = hesim_dat))
  
  # ~~ Costs ------------------
  drugcostmod <- reactive(create_StateVals(drugcost_tbl, n = n_samples(),
                                  time_reset = TRUE, hesim_data = hesim_dat))
  medcostmod <- reactive(create_StateVals(medcost_tbl, n = n_samples(),
                                 hesim_data = hesim_dat))
  costmods <- reactive(list(Drug = drugcostmod(),
                   Medical = medcostmod()))
  
  # Run the disease simulation ----------------
  
  # ~ Combining input into economic model -------------------
  
  # The above inputs are reactive() so they are updated any time the model inputs are changed.
  # The script below runs only when the input$Run_model is clicked because it is eventReactive()
  ictstm <- eventReactive(input$Run_model, {
    
    #Keeping all the hesim arguments bundled in a single expression means that the outputs will 
    # always be consistent and created at the same time
    
    #Each of the ictstm arguments in the original 'hesim example model.R' script are called here
    # as ictstm_sub, then at the end, ictstm_sub is outputted and becomes the name of the eventReactive element
    # which in this case is ictstm.
    
    ictstm_sub <- IndivCtstm$new(trans_model = transmod(),
                           utility_model = utilitymod(),
                           cost_models = costmods())
  
  # Environment objects are not functions so do not require brackets () after them - if you are not including R6 class
  #objects in your model then you will need to use brackets for all reactive data as it is now a reactive function
    
    # ~ Combining input into economic model -------------------
    
    #This shows a notification to let the user know in the front-end that the model is running, as this can take a little while
    showNotification(
      paste("Running model, this will take a few moments..."),
      id = "RunNotification",
      duration = NULL,
      type = "message"
    )
    
    # ~ Run simulation -------------
    # This runs the disease simulation, and assumed that the max patient age is 100 (after which they automatically transfer to 'Death' state)
    
  ictstm_sub$sim_disease(max_age = 100)
    
    # ~ Generate outcomes --------------
    # ~~ Survival --------------
    # Create survival curves with set time intervals
    # Time is in years, so this will measure from 0 to 30 years, with 1/12 (1 month) intervals
  ictstm_sub$sim_stateprobs(t = seq(0, n_years() , 1 / 12))
    
    # ~~ QALYS -------------
    # QALYs and costs are simulated separately from the simulation of the disease
  ictstm_sub$sim_qalys(dr = c(0,disc_QALY()))
    
    # ~~ Costs ------------
  ictstm_sub$sim_costs(dr = c(0,disc_Cost()))
    
  # The model calculations are now completed so the notification can be removed  
  removeNotification("RunNotification")
    
  # This outputs the ictstm_sub, this is like returning an object from a function. This eventReactive()
  # element that this is wrapped in is called ictstm, so this is essentially saying ictstm <- ictstm_sub
  return(ictstm_sub)
    
  })
  
  
  
   # ~ Output to shiny ---------------------------
  
  # This will output the live ictstm() results (note the brackets are added to ictstm() because it is reactive)
  # which will mean that it will update after every time the model is run. The render...() functions work like reactive()
  # functions
  output$Results_DT <- renderDataTable({
    req(input$Run_model >= 1)
    req(!is.null(ictstm()$costs_))
    
    ce_sim_ictstm <- ictstm()$summarize()
    ce_sim_ictstm <- summary(ce_sim_ictstm, labels = labs_indiv) %>%
      format()
    
    datatable(
      data =  ce_sim_ictstm,
      rownames = FALSE,
      # There are a lot of options available within datatable, such as search bars, ordering and paging options.
      # Set these to TRUE to see what they do
      options = list(
        lengthChange = FALSE,
        paging = FALSE,
        searching = FALSE,
        info = TRUE,
        ordering = FALSE,
        scrollX = FALSE,
        autoWidth = TRUE
      )
    )
  })
  
  output$Results_graph <- renderPlot({
    req(input$Run_model >= 1)
    
    autoplot(isolate(ictstm()$stateprobs_), labels = labs_indiv,
             ci = FALSE) + theme_bw()
  })
  
  output$Results_text <- renderText({
    req(input$Run_model >= 1)
    HTML("The results displayed here are based on <b>",n_patients,"</b> sampled patients with<b>",isolate(n_samples()),"</b>probabilistic samples using the input data displayed." )
  })
  
  # ~ Create report -------------------
  # This location is relative to the location of app.R
  Markdown_location <- "./"
  
  # ANSWER CODE CHANGES HERE -----------------------------
  # Model_Diagram added to the Export_params for html and pdf document
  
  output$Create_htmlreport <- downloadHandler(
      filename = "html-report_shiny.html",
      content = function(file) {
        
        # Create a notification in front-end to show this is happening
        showNotification(
          paste("Creating html report, this will take a few moments..."),
          id = "HTMLNotification",
          duration = NULL,
          type = "message"
        )
        # Make sure it closes when we exit this reactive, even if there's an error
        on.exit(removeNotification("HTMLNotification"))
        
        ce_sim_ictstm <- ictstm()$summarize()

        Export_params <- list(
          # Main results
          Stateprobs            = ictstm()$stateprobs_,
          Summarisedf           = ce_sim_ictstm,
          labs_indiv            = labs_indiv
        )

        # html document
        rmarkdown::render(
          input = file.path(Markdown_location,"html report.Rmd"),
          output_format = 'bookdown::html_document2',
          output_file = file,
          params = Export_params,
          envir = environment()
        )
      }
    )
  
  output$Create_pdfreport <- downloadHandler(
    filename = "pdf-report_shiny.pdf",
    content = function(file) {
      
      # Create a notification in front-end to show this is happening
      showNotification(
        paste("Creating pdf report, this will take a few moments..."),
        id = "PDFNotification",
        duration = NULL,
        type = "message"
      )
      # Make sure it closes when we exit this reactive, even if there's an error
      on.exit(removeNotification("PDFNotification"))
      
      ce_sim_ictstm <- ictstm()$summarize()
      
      Export_params <- list(
        # Main results
        Stateprobs            = ictstm()$stateprobs_,
        Summarisedf           = ce_sim_ictstm,
        labs_indiv            = labs_indiv
      )
      
      # pdf document
      rmarkdown::render(
        input = file.path(Markdown_location,"pdf report.Rmd"),
        output_format = 'bookdown::pdf_document2',
        output_file = file,
        params = Export_params,
        envir = environment()
      )
    }
  )
  

  #Click the model on start-up
  observe({
    if(input$Run_model == 0){
      click("Run_model")
    }
  })
   
}