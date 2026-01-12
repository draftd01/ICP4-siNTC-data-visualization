library(shiny)
library(shinyjs)
library(bslib)
library(httr)

source('home_module.R')
source('data_module.R')
source('methods_module.R')
source('analysis_module.R')
source('login_module.R')

options(shiny.session.timeout = 0)

ui <- page_fillable(
  theme = bs_theme(
    version = 5,
    navbar_bg = "#337ab7"
  ),
  
  # Add padding to account for fixed navbar
  tags$head(
    tags$style(HTML("
      .tab-content {
        padding-top: 50px !important;
      }
      .navbar-brand {
        color: white !important;
        text-decoration: none !important;
        font-size: 22px !important;
      }
      .navbar-nav .nav-item .nav-link,
      .navbar-nav .nav-item .nav-link:link,
      .navbar-nav .nav-item .nav-link:visited {
        background-color: transparent;
        color: white;
        text-decoration: none !important;
        border: none !important;
        font-size: 16px !important;
      }
      .navbar-nav .nav-item .nav-link.active,
      .navbar-nav .nav-item .nav-link.active:hover {
        background-color: white !important;
        color: black !important;
        text-decoration: none !important;
        border: none !important;
        border-bottom: none !important;
      }
      .navbar-nav .nav-item .nav-link:hover,
      .navbar-nav .nav-item .nav-link:focus {
        background-color: #d9edf7;
        color: black;
        text-decoration: none !important;
        border: none !important;
      }
      /* Add this to push login to the right */
      .navbar-nav:last-child {
        margin-left: auto;
      }
    "))
  ),
  
  navbarPage(
    title = "Secrepedia",
    id = "main_navbar",
    position = "fixed-top",
    fluid = TRUE,
    
    # Main navigation tabs
    tabPanel("Home", home_ui("home")),
    tabPanel("Data", data_ui("data")),
    tabPanel("Methods", methods_ui("methods")),
    #tabPanel("Prediction", prediction_ui("prediction")),
    #tabPanel("Analysis", analysis_ui("analysis")),
    
    # Login dropdown on the right
    #login_ui("login"),
    #removing annoying styling from these elements
    tags$script(HTML("
      $(document).ready(function() {
        console.log($('.navbar-header'))
        console.log($('#main_navbar'))
        $('.navbar-header').removeClass('navbar-header');
        $('#main_navbar').css('flex-direction', 'unset');
      });
    "))
  )
)

#THIS FUNCTION INITIALIZES THE SERVER
server <- function(input, output, session) {
  # First initialize login module
  experiments <- get_experiments()#reactiveValues(names=get_experiments())# get_experiments()#login_server("login", session)
  #browser()
  #print("INITIALIZING SERVER")
  # Then call the data module server
  data_server("data", session, experiments)
  #prediction_server("prediction")
  
}

# Run the Shiny App
shinyApp(ui, server)