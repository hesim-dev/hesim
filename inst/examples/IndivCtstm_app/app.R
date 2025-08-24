#The app.R script is where the app can be called from
# It can also be called from hesim::shinyIndivCtstm()
# See https://shiny.rstudio.com/gallery/ for examples of shiny app layouts and functionalties

library("shiny")
library("hesim")
options(encoding = "UTF-8")

# Setwd ---------------
# setwd if running script from this file. Set wd() as the hesim project location
# setwd()

App_location <- "./inst/examples/IndivCtstm_app/"

runApp(appDir = file.path(getwd(),App_location),     # This looks for the ui.R, server.R and www folder within the wd. Do not rename these.
       launch.browser = TRUE,                    # This line will open the application in the user's default browser
       quiet = FALSE,
       display.mode = "normal")
