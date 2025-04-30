library(shiny)
library(bslib)

login_ui <- function(id) {
  ns <- NS(id)
  
  # Dropdown menu items for login using nav_menu
  nav_menu(
    title = "Login",
    align = "right",
    nav_panel(
      value = "login_panel",
      tags$div(
        class = "dropdown-menu dropdown-menu-end p-3 show",
        style = "min-width: 300px; position: static;",
        draggable = "false",  # Prevent dragging at container level
        
        # Welcome message container that will be updated
        tags$div(
          class = "text-center mb-3",
          uiOutput(ns("welcome_message"))
        ),
        
        tags$div(
          style = "position: relative;",
          draggable = "false",  # Prevent dragging at inner container level
          tags$input(
            id = ns("username"),
            type = "text",
            class = "form-control",
            placeholder = "Enter username",
            style = "user-select: text !important; -webkit-user-select: text !important;",
            onkeydown = "event.stopPropagation();",
            onmousedown = "event.stopPropagation();",
            onmousemove = "event.stopPropagation();",
            onclick = "event.stopPropagation();",
            ondragstart = "return false;",  # Prevent drag start
            draggable = "false"  # Prevent dragging at input level
          ),
          tags$input(
            id = ns("password"),
            type = "password",
            class = "form-control mt-3",
            placeholder = "Enter password",
            style = "user-select: text !important; -webkit-user-select: text !important;",
            onkeydown = "event.stopPropagation();",
            onmousedown = "event.stopPropagation();",
            onmousemove = "event.stopPropagation();",
            onclick = "event.stopPropagation();",
            ondragstart = "return false;",  # Prevent drag start
            draggable = "false"  # Prevent dragging at input level
          )
        ),
        div(
          class = "d-grid gap-2 mt-3",
          actionButton(ns("login_button"), "Login", class = "btn-primary"),
          div(
            class = "text-danger",
            textOutput(ns("login_message"))
          )
        )
      )
    )
  )
}

login_server <- function(id, parent_session) {
  moduleServer(id, function(input, output, session) {

    user_info <- reactiveValues(is_authenticated = FALSE, username = NULL, token = NULL, experiments = NULL)
    #url_base = "https://api.secrepedia.org/secrepediadb/"
    url_base = "https://dev-api.secrepedia.org/secrepediadb/"
    
    # Render welcome message
    output$welcome_message <- renderUI({
      if (user_info$is_authenticated && !is.null(user_info$username)) {
        tags$p(
          style = "font-size: 0.9em; margin: 0;",
          paste("Logged in as", user_info$username)
        )
      } else {
        tagList(
          #tags$h6("Welcome to Secrepedia"),
          tags$p("Please login to access all features")
        )
      }
    })
                                          
    observeEvent(input$login_button, {
    
    # Define the API URL for Django authentication      
    user_info$username <- input$username
    user_info$password <- input$password

    if (is.null(user_info$username) || is.null(user_info$password)){
      output$login_message <- renderText("Guest user, only public experiments are available.")
      user_info$is_authenticated <- TRUE
      user_info$access_token <- NULL
      user_info$refresh_token <- NULL
      # Add shorter delay and switch to data tab for guest login
      invalidateLater(500, session)
      observe({
        shinyjs::runjs("document.querySelector('a[data-value=\"Data\"]').click();")
      })
      
    }
    else{
      # Prepare the payload to be sent to Django
      payload <- list(
          username = user_info$username,
          password = user_info$password)
      # Send POST request to Django backend for login
      response <- httr::POST(paste0(url_base, "token/"), body = payload, encode = "json")

      if (response$status_code == 200) {
          # Parse response and store session info
          tokens <- httr::content(response)
          user_info$is_authenticated <- TRUE
          user_info$access_token <- tokens$access
          user_info$refresh_token <- tokens$refresh
          output$login_message <- renderText("Login successful!")
          # Add shorter delay and switch to data tab after successful login
          invalidateLater(500, session)
          observe({
            shinyjs::runjs("document.querySelector('a[data-value=\"Data\"]').click();")
          })
        }
        else{
          # Feedback on login failure and switch to data tab
          user_info$is_authenticated <- FALSE
          output$login_message <- renderText("Login failed. Only public experiments are available.")
          user_info$access_token <- NULL
          user_info$refresh_token <- NULL
          invalidateLater(500, session)
          observe({
            shinyjs::runjs("document.querySelector('a[data-value=\"Data\"]').click();")
          })
    } 
    }
      
    })
    
    observe({
    # Fetch experiments
    if (is.null(user_info$access_token)){
      # send request without token
      exp_response <- httr::GET(
        paste0(url_base, "experiments/"), config = list())
     
    }else{
       exp_response <- httr::GET(paste0(url_base, "experiments/"),
                          httr::add_headers(Authorization = sprintf("Bearer %s", user_info$access_token)))
    }

    if (exp_response$status_code == 200){
        user_info$experiments <- httr::content(exp_response, "parsed", "application/json")
        output$exp_message <- renderText("Experiments fetched successfully! Please go to the data tab to view experiments.")
    } 
    else if (exp_response$status_code == 401){
        # try refresh token
        refresh_response <- httr::POST(paste0(url_base, "token/refresh/"),
                                       body = list(refresh = user_info$refresh_token),
                                       encode = "json")
        if (refresh_response$status_code == 200){
            tokens <- httr::content(refresh_response)
            user_info$access_token <- tokens$access
            # try again
            exp_response <- httr::GET(paste0(url_base, "experiments/"),
                          httr::add_headers(Authorization = sprintf("Bearer %s", user_info$access_token)))
            if (exp_response$status_code == 200){
                user_info$experiments <- httr::content(exp_response, "parsed", "application/json")
                output$exp_message <- renderText("Experiments fetched successfully! Please go to the data tab to view experiments.")
            }
            else{
                output$exp_message <- renderText("Token fail to refresh. \nPlease login again.")
            }
        }
        else{
            output$exp_message <- renderText("Token fail to refresh. \n Please login again.")
        }
    } else {
        print(paste("Unexpected status code:", exp_response$status_code))
    }
    })
    
    # Return reactive values for use in main app
    return(user_info)
  })
}
