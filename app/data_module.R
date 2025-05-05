library(bslib)
source('plot_module.R')

data_ui <- function(id) {
  ns <- NS(id)
  tabsetPanel(
    id = ns("data_tabs"),
    tabPanel(
      "Overview",
      layout_column_wrap(
        width = 1,
        heights_equal = "row",
        card(
          card_header("Available Datasets"),
          uiOutput(ns("experiment_links"))
        ),
        uiOutput(ns("exp_cards_container"))
      )
    )
    # Additional tabs will be inserted here
  )
}

data_server <- function(id, main_session, user_info) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    output$experiment_links <- renderUI({
      # Create ordered experiment list once to ensure consistency
      print("experiments")
      #print(user_info$experiments)
      ordered_experiments <- user_info$experiments[order(sapply(user_info$experiments, function(exp) exp$id))]
      print(dim(search_results()))
      print(length(ordered_experiments))
      div(
        style = "display: flex; align-items: start; width: 100%;",
        # Left column with experiment links
        div(
          style = "min-width: fit-content;",
          tags$div(style = "height: 25px;"),
          tags$ul(
            style = "margin: 0; padding: 0; list-style-type: none;",
            lapply(ordered_experiments, function(exp) {
              tags$li(
                style = "height: 20px; display: flex; align-items: center; 
                         padding-right: 20px; margin: 2px 0;",
                  actionLink(ns(paste0("exp_link_", exp$id)), exp$title)
                )
            })
          )
        ),
        # Middle column with gene grid (when present)
        if (!is.null(search_results()) && dim(search_results())[1] == length(ordered_experiments)) {
          div(
            style = "position: sticky; left: 0;",
            createGeneGrid(ordered_experiments, search_results())
          )
        },
        # Right column with search panel
        div(
          style = "min-width: 200px; padding-left: 20px; border-left: 1px solid #ddd; 
                   margin-left: auto;",
          h5("Quick molecule search"),
          selectInput(ns("id_type"), "Select ID Type",
                     choices = c(
                       "Gene Name" = "gene_name",
                       "Ensembl Gene ID" = "ensembl_gene_id",
                       "Ensembl Peptide ID" = "ensembl_peptide_id",
                       "Uniprot ID" = "uniprot"
                     ),
                     selected = "gene_name"
          ),
          textInput(ns("gene_search"), "IDs to search", 
                   placeholder = "Enter IDs (comma-separated)"),
          div(
            style = "display: flex; gap: 10px;",
            actionButton(ns("search_btn"), "Search", icon = icon("magnifying-glass")),
            actionButton(ns("clear_btn"), "Clear", icon = icon("xmark"))
          )
        )
      )
    })
    
    # Reactive values to store active experiment cards and selected experiment
    active_cards <- reactiveValues(cards = list())
    created_tabs <- reactiveValues(tabs = list())
    #selected_exp <- reactiveVal(NULL)  # Store the selected experiment
    
    # Reactive elements and helper functions to show presence of genes in experiments
    search_results <- reactiveVal(NULL)
    
    observeEvent(input$search_btn, {
      req(input$gene_search, input$id_type)
      genes <- trimws(strsplit(input$gene_search, ",")[[1]])
      
      # Use the same ordered experiments
      ordered_experiments <- user_info$experiments[order(sapply(user_info$experiments, function(exp) exp$id))]
      exp_ids_to_search <- as.character(sapply(ordered_experiments, function(exp) exp$id))
      
      #print(paste("exp_ids_to_search:", exp_ids_to_search))
      print(paste("user_token NULL:", is.null(user_info$access_token)))
      
      # Assuming presence_matrix is a data frame where:
      # - row names are experiment IDs
      # - column names are gene IDs
      # - values are TRUE/FALSE for presence/absence
      results <- search_molecules(user_info, input$id_type, genes, exp_ids_to_search)
      print("search_molecules results status:")
      print(results$status)
      print(results$msg)
      presence_matrix <- results$out
    
      print(presence_matrix)
      # display warning message if colnames(presence_matrix) is different from genes
      if (is.null(presence_matrix) || length(setdiff(genes, colnames(presence_matrix))) != 0) {
        showModal(modalDialog(
          title = "",
          "Some genes were not found in the database. Please check id type or spelling."
        ))
      }
      if (is.null(presence_matrix)) {
        search_results(NULL)
      }
      else{
      display_presence <- matrix(FALSE, nrow = length(exp_ids_to_search), ncol = dim(presence_matrix)[2])
      rownames(display_presence) <- exp_ids_to_search
      colnames(display_presence) <- colnames(presence_matrix)
      display_presence[match(rownames(presence_matrix), exp_ids_to_search), ] <- presence_matrix
      print(display_presence)
      search_results(display_presence)
      }
    })
    
    # Add clear button observer
    observeEvent(input$clear_btn, {
      search_results(NULL)  # Reset the search results
      updateTextInput(session, "gene_search", value = "")  # Clear the input field
    })
    
    # Helper function for the gene grid
    createGeneGrid <- function(experiments, gene_results) {
      # Error handling wrapper
      tryCatch({
        tags$div(
          style = "display: inline-block;",
          # Gene names as column headers
          tags$div(
            style = "display: grid; grid-auto-columns: auto; grid-auto-flow: column; gap: 40px; margin-bottom: 10px;",
            tags$div(style = "width: 40px;"),
            lapply(colnames(gene_results), function(gene) {
              tags$div(
                style = "display: flex; justify-content: center; align-items: center; 
                         min-width: 60px; padding: 0 10px; height: 20px; font-weight: bold;",
                gene
              )
            })
          ),
          # Grid of dots
          tags$div(
            lapply(experiments, function(exp) {
              tags$div(
                style = "display: grid; grid-auto-columns: auto; grid-auto-flow: column; gap: 40px; 
                         height: 20px; align-items: center; margin: 2px 0;",
                tags$div(style = "width: 40px;"),
                lapply(colnames(gene_results), function(gene) {
                  is_present <- gene_results[as.character(exp$id), gene]
                  tags$div(
                    style = "display: flex; justify-content: center; align-items: center; 
                             min-width: 60px; padding: 0 10px; height: 20px;",
                    tags$div(
                      style = sprintf("width: 12px; height: 12px; border-radius: 50%%; 
                                     background-color: %s; opacity: 0.8;",
                                     if(is_present) "#ff1493" else "#f5f5f5")
                    )
                  )
                })
              )
            })
          )
        )
      },
      error = function(e) {
        # On error, return just the experiment links
        tags$div(
          style = "color: red; margin-bottom: 10px;",
          paste("Error:", e$message)
        )
      })
    }

  
    
    # Utility function to generate a card UI dynamically
    generate_exp_card <- function(exp_id, exp_title, user_info) {
      div(
        class = "mb-3",  # adds margin-bottom of size 3 (Bootstrap spacing)
        card(
          id = ns(paste0("card_", exp_id)), 
          card_header(exp_title),
          fluidRow(
            column(6, show_exp_details(user_info, exp_id)),
            #column(4, actionButton(ns(paste0("get_tab_", exp_id)), "Get Analysis Results")),
            #column(2, actionButton(ns(paste0("remove_", exp_id)), "Hide"))
            column(4, div(style = "text-align: right;",
              actionButton(ns(paste0("get_tab_", exp_id)), "Get Analysis Results")
            )),
            column(2, div(style = "text-align: right;",
              actionButton(ns(paste0("remove_", exp_id)), "Hide")
            ))
          )
        )
      )
    }
    
    # Render UI for cards
    output$exp_cards_container <- renderUI({
      do.call(tagList, active_cards$cards)
    })
    
     # Wrap the experiment link and remove button setup in observe()
    observe({
      lapply(user_info$experiments, function(exp) {
        exp_link_id <- paste0("exp_link_", exp$id)
        
        observeEvent(input[[exp_link_id]], {
          if (!as.character(exp$id)%in% names(active_cards$cards)) {
            active_cards$cards[[as.character(exp$id)]] <- generate_exp_card(exp$id, exp$title, user_info)
            #selected_exp(exp$id)
          }
        })
        
        # Set up remove button observer for this experiment
        observeEvent(input[[paste0("remove_", exp$id)]], {
          active_cards$cards[[as.character(exp$id)]] <- NULL
        })
      })
    })

    # Helper function to insert or switch to a tab
    insert_or_switch_tab <- function(exp_id, exp_name, molecule_type) {
      tab_id <- paste0("tab_", gsub(" ", "_", exp_id))
      close_button_id <- paste0("close_", tab_id)
      plot_id <- paste0("plot_", tab_id) 
      print(paste("Current input$data_tabs value:", input$data_tabs))  # Debug print
      if (!tab_id %in% names(created_tabs$tabs)) {
        print(paste("Creating new tab", tab_id))
        new_tab <- tabPanel(
          title = paste("Experiment", exp_id),
          value = tab_id,
          fluidRow(
            class = "mt-2 mb-4",  # mt-2 reduces top margin, mb-4 increases bottom margin
            column(9, h4(exp_name, class = "mb-0")),  # mb-0 removes margin below h3
            column(3, actionButton(ns(close_button_id), "Close Tab"))
          ),
          plot_ui(ns(plot_id))
        )
        
        insertTab(
          session = session,
          inputId = "data_tabs",
          tab = new_tab,
          target = "Overview",
          position = "after",
          select = TRUE
        )
        
        # Store the module instance
        created_tabs$tabs[[tab_id]] <- list(
          name = exp_name,
          plot_id = plot_id,
          module = plot_server(plot_id, user_info, exp_id, molecule_type)
        )
        
        # Set up remove button observer for this tab
        observeEvent(input[[close_button_id]], {
          print(paste("Closing tab:", tab_id))
          
          # Get the plot module's local_selected and trigger same behavior as clear_genes
          if (!is.null(created_tabs$tabs[[tab_id]]$module)) {
            # Get the module instance
            module <- created_tabs$tabs[[tab_id]]$module
            
            # Trigger same behavior as clear_genes
            module$local_selected$should_render_heatmap <- FALSE
            module$local_selected$genes <- character(0)
            module$local_selected$click <- character(0)
            module$local_selected$typed <- character(0)
            
            # Then do the regular cleanup
            module$cleanup()
          }
          
          removeTab(session = session, inputId = "data_tabs", target = tab_id)
          created_tabs$tabs[[tab_id]] <- NULL
          updateTabsetPanel(session, "data_tabs", selected = "Overview")
        }, ignoreInit = TRUE)
        
      } else {
        print(paste("Switching to existing tab:", tab_id))
        updateTabsetPanel(session, "data_tabs", selected = tab_id)
      }
    }

    # Observer for the Get Analysis Results button, open or switch to the tab
    observe({
      lapply(user_info$experiments, function(exp) {
        tab_id <- paste0("get_tab_", exp$id)
        
        observeEvent(input[[tab_id]], {
          insert_or_switch_tab(exp$id, exp$title, exp$molecule_type)
        }, ignoreInit = TRUE, ignoreNULL = TRUE)
      })
    })
  })
}
