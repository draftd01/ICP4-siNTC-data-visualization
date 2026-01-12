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
    ),
    # Additional tabs will be inserted here
    
    #-----STYLES-----
    tags$head(
      tags$style(HTML("
         .legend-item{
            display: flex;
            align-items: center;
            gap: 5px;
            font-size: 0.8rem;
         }  
         
         .dot{
            width: 12px;
            height: 12px;
            border-radius: 50%;
         }
         
         .link{
            height: 20px;
            display: flex;
            align-items: center; 
            padding-right: 20px;
            margin: 2px 0;
         }
      "))
    )
    #-----END STYLES-----
  )
}

data_server <- function(id, main_session, experiments) { #dir names of experiments in "experiments" folder, format <exp>_comparisons
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    #-----INITILIZE REACTIVE VALUES-----
    
    experiment_data <- reactiveValues(
      name = NULL,
      #txt_files = NULL, #all txt files selected experiment folder with ext (regardless of filter)
      selected_files = NULL, #txt files without ext (could be filtered or not)
      all_files = NULL, #list of all experiment txt files (without ext)
      all_files_full_path = NULL
    )
    # Reactive values to store active experiment cards and selected experiment
    active_card <- reactiveVal(NULL)
    created_tabs <- reactiveValues(tabs = list())
    links <- reactiveValues(observers = list())
    cards <- reactiveValues(observers = list())
    previously_selected <- reactiveValues(exps = list())
    search_params <- reactiveValues(p_val = 0.05, lfc_cutoff = 0, genes=NULL)
    search_results <- reactiveVal(NULL)
    gene_data <- reactiveVal(NULL)
    
    
    #-----END INITIALIZE REACTIVE VALUES-----
    
    #-----PAGE RENDER SECTION-----
    
    output$experiment_links <- renderUI({
      # Create ordered experiment list once to ensure consistency
      
      #experiments <- experiment_names$names
      #browser()
      # <- user_info$experiments
      #base_txt_files <- NULL

      div(
        style = "display: flex; align-items: start; width: 100%;",
        
        # Left column with experiment links
        if(!is.null(experiment_data$all_files)){ 
          div(
            style = "min-width: fit-content; display:flex; overflow-x:auto;",
            tags$ul( #list of all exp
              style = "margin: 0; padding: 0; list-style-type: none;",
              lapply(experiments, function(exp_name) {
                exp_label <- strsplit(exp_name, "_")[[1]][1] #exp name without "_comparisons"
                tags$li( #exp list item
                  class="pt-3",
                  div(
                    class="d-flex",
                    div(
                      style="min-width: fit-content;",
                      HTML(paste("<b><u>",exp_label,"</u></b>")),
                      tags$ul(#list of all files under each exp
                        lapply(list_txt_files(exp_name),function(file_name){
                              tags$li(class="link", #file list item
                                      actionLink(ns(file_name),file_name))
                          }
                        ),
                        if(!is.null(search_results())){
                          tags$li( 
                            class = "d-block",
                            "ChipSeq",
                          )
                        }
                      )
                    ),
                    div(
                      if(!is.null(search_results())){
                        createGeneGrid(search_results(),exp_label=="siNTC")
                      }
                    )
                    
                  )
                )
              }),
              # Legend section
              if(!is.null(search_results())){
                div(class="d-flex",
                  div(class="w-100"),#filler div
                  div(
                    style = "display: flex; justify-content: center; column-gap: 20px; margin-top: 15px; flex-wrap:wrap; min-width:225px; margin-right:15px; margin-left:auto;",
                    div(
                      class = "legend-item",
                      div( class = "dot",
                           style = "background-color: dodgerblue;"),
                      "Downregulated"
                    ),
                    div(
                      class = "legend-item",
                      div( class = "dot",
                           style = "background-color: darkorange;"),
                      "Upregulated"
                    ),
                    div(
                      class = "legend-item",
                      div( class = "dot",
                           style = "background-color: lightgray;"),
                      "Not significant"
                    ),
                    div(
                      class = "legend-item",
                      div( class = "dot",
                           style = "background-color: slategray;"),
                      "NA"
                    ),
                    div(
                      class = "legend-item",
                      div(class = "dot",
                          style = "background-color: limegreen;"),
                      "Found"
                    ),
                    div(
                      class = "legend-item",
                      div(class = "dot",
                          style = "background-color: red;"),
                      "Not found"
                    )
                  )
                )
              }
            )
          )
        },
        div(
          style = "min-width: 200px; padding-left: 20px; border-left: 1px solid #ddd; 
                 margin-left: auto;",
          h5("Quick molecule search"),
          selectInput(ns("id_type"), "Select ID Type",
                      choices = c(
                        "Gene Name" = "gene_name", #first is text in drop down, second is actual value submitted in form
                        "Ensembl Gene ID" = "ensembl_gene_id"
                        #"Ensembl Peptide ID" = "ensembl_peptide_id",
                        #"Uniprot ID" = "uniprot"
                      ),
                      selected = "gene_name"
          ),
          textInput(ns("gene_search"), "Values to search", 
                    placeholder = "Enter values (comma-separated)", value=search_params$genes),
          numericInput(ns("p_val_gene_grid"),"p-value cutoff", search_params$p_val),
          numericInput(ns("lfc_cutoff_gene_grid"), "Log fold cutoff", search_params$lfc_cutoff),
          div(
            style = "display: flex; gap: 10px;",
            actionButton(ns("search_btn"), "Search", icon = icon("magnifying-glass")),
            actionButton(ns("clear_btn"), "Clear", icon = icon("xmark"))
          )
        )
      )
    })
    
    # Render UI for cards
    output$exp_cards_container <- renderUI({
      #print("rendering active card")
      req(active_card())
      active_card()
    })
    
    #-----END PAGE RENDER SECTION-----

    # -----REACTIVE ELEMENTS/HELPER FUNCTIONS TO FILTER EXPERIMENTS BY PRESENCE OF SPECIFED GENES-----
    
    observeEvent(input$search_btn, {
      req(input$gene_search, input$id_type) #get input and its type (gene name, id, etc)
      
      genes <- trimws(strsplit(input$gene_search, ",")[[1]]) 
      search_params$p_val <- input$p_val_gene_grid
      search_params$lfc_cutoff <- input$lfc_cutoff_gene_grid
      search_params$genes <- input$gene_search
      
      # Apresence_matrix is df where:
      # - row names are experiment IDs
      # - column names are gene IDs
      # - values are TRUE/FALSE for presence/absence
      full_folder_name <- paste0( c("experiments/", experiment_data$name, "_comparisons"), collapse="" )
      presence_df <- search_molecules(input$id_type, genes, experiment_data$all_files_full_path, input$p_val_gene_grid, input$lfc_cutoff_gene_grid) #returns df with rownames as genes and colnames as experiments; values either TRUE/FALSE for gene presence
      # display warning message if colnames(presence_matrix) is different from genes
      # if (is.null(presence_df) || length(rownames(presence_df)) < length(genes)) {
      #   showModal(modalDialog(
      #     title = "",
      #     "Some genes were either filtered out during quality control or not found in the database. Please check id type or spelling."
      #   ))
      # }
      
      if (is.null(presence_df)) {
        search_results(NULL)
      }
      else{
        display_presence <- as.matrix(presence_df)
        search_results(display_presence)
      }
      
    })
    
    #-----END FILTER FUNCTIONS-----
    
    # Add clear button observer
    observeEvent(input$clear_btn, {
      search_results(NULL)  # Reset the search results
      experiment_data$selected_files <- gsub("\\.txt$","",experiment_data$txt_files)
      updateTextInput(session, "gene_search", value = "")  # Clear the input field
    })
    
    # Helper function for the gene grid
    createGeneGrid <- function(gene_results, siNTC) { #gene_results: rownames = genes, colnames = experiments
      # Error handling wrapper
      exp_files <- colnames(gene_results)[!xor(siNTC,grepl("siNTC",colnames(gene_results)))]
      #print(paste("main exp name",exp_name))
      #print(colnames(gene_results))
      #print(rownames(gene_results))
      tryCatch({
        div(
          style = "display: inline-block;",
          # Gene names as column headers
          #div(style = "height: 25px;"),
          div( #gene grid column headers
            style = "display: grid; grid-auto-columns: auto; grid-auto-flow: column; gap: 40px; margin-bottom: 10px;",
            lapply(rownames(gene_results), function(gene) {
              div(
                style = "display: flex; justify-content: center; align-items: center; 
                         min-width: 60px; padding: 0 10px; height: 20px; font-weight: bold;",
                gene
              )
            })
          ),
          # Grid of dots
          div(
            lapply(exp_files, function(file) { #file = file_name (no ext, just as in presence matrix)
              #print(file) #ensuring files are in same order as list
              div( #table row
                style = "display: grid; grid-auto-columns: auto; grid-auto-flow: column; gap: 40px; 
                         height: 20px; align-items: center; margin: 2px 0;",
                lapply(rownames(gene_results), function(gene) {
                  lfc <- gene_results[gene, file]
                  div( #table cell
                    style = "display: flex; justify-content: center; align-items: center; 
                             min-width: 60px; padding: 0 10px; height: 20px;",
                    div(
                      style = sprintf("width: 12px; height: 12px; border-radius: 50%%; 
                                     background-color: %s; opacity: 0.8;",
                                     switch(as.character(lfc), "-1" = "dodgerblue", "0" = "lightgray", "1" = "darkorange", "-2" = "red", "2" = "limegreen", "slategray")) #,downregulated, ns, or upregulated, not found, or NA (no p value)
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
        div(
          style = "color: red; margin-bottom: 10px;",
          paste("Error:", e$message)
        )
      })
    }
    
    # Utility function to generate a card UI dynamically
    generate_exp_card <- function(file_name) { #file_name no ext.
      print(paste("card:", file_name))
      #base_file_name <- gsub("\\.txt$", "", file_name)
      div(
        class = "mb-3",  # adds margin-bottom of size 3 (Bootstrap spacing)
        card(
          id = ns(paste0("card_", file_name)), 
          card_header(file_name),
          fluidRow(
            column(6, div()),#show_exp_details(user_info, exp_id)),
            column(4, div(style = "text-align: right;",
              actionButton(ns(paste0("get_tab_", file_name)), "Get Analysis Results")
            )),
            column(2, div(style = "text-align: right;",
              actionButton(ns(paste0("remove_", file_name)), "Hide")
            ))
          )
        )
      )
    }
    
     # Experiment onclick event to bring up txt file list under this experiment directory
    # observe({
    #   #experiments <- experiment_names$names
    #   lapply(experiments, function(exp){ #
    # 
    #     exp_title <- strsplit(exp, "_")[[1]][1] #exp folder names are of format exp_comparisons
    # 
    #     observeEvent(input[[exp_title]], {
    # 
    #       experiment_data$name <- exp_title #without "_comparison"
    # 
    #       # Find corresponding full folder name
    #       txt_files_path <- file.path("experiments", exp)
    # 
    #       if (dir.exists(txt_files_path)) {
    #         txt_files <- list.files(path = txt_files_path, pattern = "*\\.txt$", full.names = FALSE)
    #         txt_files <- txt_files[!grepl("^(GSEA|logCPM)", txt_files)] #only get non-GSEA/logCPM data (only DEseq)
    #         txt_files <- txt_files[!grepl("Shrink", txt_files)]
    #         experiment_data$txt_files <- txt_files
    #         base_txt_files <- gsub("\\.txt$","",txt_files)
    #         experiment_data$selected_files <- base_txt_files
    #       } else {
    #         experiment_data$txt_files <- NULL
    #         experiment_data$selected_files <- NULL
    #       }
    #     })
    #   })
    # })
    # observe({
    #   print("INSIDE THE EVENT OBSERVER")
    #   exp_title <- strsplit(exp, "_")[[1]][1] #exp folder names are of format exp_comparisons
    #   experiment_data$name <- exp_title #without "_comparison"
    #   
    #   # Find corresponding full folder name
    #   txt_files_path <- file.path("experiments", exp)
    #   
    #   if (dir.exists(txt_files_path)) {
    #     txt_files <- list.files(path = txt_files_path, pattern = "*\\.txt$", full.names = FALSE) 
    #     txt_files <- txt_files[!grepl("^(GSEA|logCPM)", txt_files)] #only get non-GSEA/logCPM data (only DEseq)
    #     txt_files <- txt_files[!grepl("shrink", txt_files)]
    #     experiment_data$txt_files <- txt_files #files with ext
    #     #set link observers
    #     links$file_observers <- lapply(txt_files, function(file_name) {
    #       base_file_name <- gsub("\\.txt$","",file_name)
    #       observeEvent(input[[base_file_name]], {
    #         print(paste("active card for", base_file_name, collapse=" "))
    #         active_card(generate_exp_card(file_name))
    #       })
    #     })
    #     base_txt_files <- gsub("\\.txt$","",txt_files)
    #     experiment_data$selected_files <- base_txt_files #files w/o ext
    #     
    #   } else {
    #     experiment_data$txt_files <- NULL
    #     experiment_data$selected_files <- NULL
    #     #print("Directory not found or no txt files")
    #   }
    # })


    # Helper function to insert or switch to a tab, also set focus of selected data file
    insert_or_switch_tab <- function(full_file_path){
      base_file_name <- gsub("\\.txt$", "", basename(full_file_path))
      base_file_name <- gsub("-", "_", base_file_name) #used for CARD HEADER, TAB ID, CLOSE BUTTON ID, AND PLOT ID
      tab_id <- paste0("tab_", base_file_name)
      close_button_id <- paste0("close_", tab_id)
      plot_id <- paste0("plot_", tab_id) 

      if (!tab_id %in% names(created_tabs$tabs)) {
        new_tab <- tabPanel(
          title = paste("Experiment", base_file_name),
          value = tab_id,
          fluidRow(
            class = "mt-2 mb-4",  # mt-2 reduces top margin, mb-4 increases bottom margin
            column(9, h4(base_file_name, class = "mb-0")),  # mb-0 removes margin below h3
            column(3, actionButton(ns(close_button_id), "Close Tab"))
          ),
          plot_ui(ns(plot_id)) #plots log2FC data from file
        )
        
        insertTab(
          session = session,
          inputId = "data_tabs",
          tab = new_tab,
          target = "Overview",
          position = "after",
          select = TRUE
        )
        
        #browser()
        # Store the module instance
        created_tabs$tabs[[tab_id]] <- list( #tab_<base_file_name>
          name = base_file_name,
          plot_id = plot_id,
          module = plot_server(plot_id, experiments, full_file_path, gene_data()) #experiments was originally user_info
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
        
      } else { # switching to new tab
        #print(paste("Switching to existing tab:", tab_id))
        updateTabsetPanel(session, "data_tabs", selected = tab_id)
      }
    } # ----END insert_or_switch_tab----
    
    list_txt_files <- function(exp_name,full_path=FALSE){ #(no path), returns filtered list of files in directory (no ext.)
      # Find corresponding full folder name
      #print("printing exp_name")
      #print(exp_name)
      txt_files_path <- file.path("experiments", exp_name)
      
      if (dir.exists(txt_files_path)) {
        txt_files <- list.files(path = txt_files_path, pattern = "*\\.txt$", full.names = full_path) 
        txt_files <- txt_files[!grepl("(GSEA|logCPM)", txt_files)] #only get non-GSEA/logCPM data (only DEseq)
        txt_files <- txt_files[!grepl("Shrink", txt_files)]
        if(!full_path)
          txt_files <- gsub("\\.txt$","",txt_files)
        return(txt_files)
      }
    }

    # Observer for the Get Analysis Results button, open or switch to the tab
    observe({
      lapply(experiment_data$all_files_full_path, function(full_file_path) { # w/ ext.
        base_file_name <- gsub("\\.txt$","", basename(full_file_path))
        tab_id <- paste0("get_tab_", base_file_name)
        observeEvent(input[[tab_id]], {
          insert_or_switch_tab(full_file_path)
        }, ignoreInit = TRUE, ignoreNULL = TRUE)
      })
    })
    
    # Implement back button
    observeEvent(input$back_to_experiments, {
      experiment_data$name <- NULL
      #experiment_data$txt_files <- NULL
      experiment_data$selected_files <- NULL
      search_results(NULL)
      active_card(NULL)
      #links$file_observers <- list()
    })
    
    observe({
      print("observing event")
      if(is.null(experiment_data$all_files)){
        #print("printing all experiments")
        #print(experiments)
        experiment_data$all_files <- unlist(lapply(experiments,list_txt_files))
        experiment_data$all_files_full_path <- unlist(lapply(experiments,list_txt_files,full_path=TRUE))
        gene_info <- c()
        for(file in experiment_data$all_files_full_path){
          df <- read.table(file, header=TRUE)
          gene_ids <- df[["gene_id"]]
          unlabeled_genes <- setdiff(gene_ids,names(gene_info))
          filtered_df <- df[df[["gene_id"]] %in% unlabeled_genes,]
          gene_info[filtered_df[["gene_id"]]] <- filtered_df[["gene"]]
        }
        gene_data(gene_info)
        #browser()
        #print(experiment_data$all_files_full_path)
        req(experiment_data$all_files)
      #if(!(experiment_data$name %in% previously_selected$exps)){#prevent observer accumulation by checking if experiment already previously selected before assinging new observers
        #previously_selected$exps <- append(previously_selected$exps, experiment_data$name)
        lapply(experiment_data$all_files, function(file_name) { #no ext
          observeEvent(input[[file_name]], {
            #req(experiment_data$name)      # double-check we did not leave
            active_card(NULL)
            active_card(generate_exp_card(file_name)) # no ext
          }, ignoreInit = TRUE)
          
          observeEvent(input[[paste0("remove_", file_name)]], { #no ext
            active_card(NULL)
          }, ignoreInit = TRUE)
        })
      }
    })
  })
}
