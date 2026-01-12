library(base)
library(shiny)
library(plotly)
library(DT)

source('utils.R')

# Combined Plot Module: plot_module.R

plot_ui <- function(id){
  #dat_dir = '/gpfs/data/proteomics/projects/fat_secretome/pipeline/data/'
  #dat_dir = '/srv/shiny-server/data/'
  
  #ID is full file path
  print(paste0("plot_ui id: ", id))
  ns <- NS(id)
  
  # Global styling applied at the top level
  tagList(
    tags$style(HTML("
      /* Ensure cards expand with content */
      .card {
        height: auto !important;
        overflow: visible !important;
      }
      
      /* Set a flexible height for the plot column */
      .plot-column {
        height: auto !important;
      }
      
      /* Additional style for plotly output to help with resizing */
      .plotly-output {
        height: auto;
      }

      .modebar-container {
        position: relative;
        top: -25px !important; 
      }

    ")),
    tags$head(
        tags$style(HTML("
            .shiny-notification {
                position: fixed;
                top: 50%;
                left: 50%;
                transform: translate(-50%, -50%);
                width: 300px;
            }
        "))
    ),
    
    card(
      card_header(h5('Differential Analysis')),
      fluidRow(
        column(4,
             selectInput(ns("comparison"), "Select Comparison:", 
                         choices = NULL),
             sliderInput(ns("p_cutoff"), "-log10(P) Cutoff:", 
                         min = 0, max = 10, value = 2),
             sliderInput(ns("lfc_cutoff"), "LFC Cutoff:", 
                         min = 0, max = 10, value = 2),
             selectInput(ns("sort_by"), "Sort Table By:", 
                         choices = c("Positive Log2FC", "Negative Log2FC", "P Value", "Adjusted P Value")),
             numericInput(ns("table_cutoff"), "Cutoff for Sorted: ",
                          value = 0),
             numericInput(ns("n_feature"), "Number of Features:", 
                          min = 0, max = 500, step = 1, value = 20)
        ),
        column(8, class = "plot-column",
           card(
             style = "width: 100%; padding-top: 20px;",
             plotlyOutput(ns("volcano_plot"), 
                         height = "calc(100vw * 0.35)",
                         width = "100%"
             )
           )
        ),
      ),
      fluidRow(
        column(8, 
               DTOutput(ns("top_table"))  # Top table
        )
      )
    ),
    
    card(
      card_header(h5("Normalized Mean Log Intensity")),
      fluidRow(
        column(4, 
               selectInput(ns("id_type"), "id Type:", choices = c(
                 "Gene Name" = "gene_name", #first is text in drop down, second is actual value submitted in form
                 "Ensembl Gene ID" = "ensembl_gene_id"
               )),
               textAreaInput(ns("genes"), "Selected Genes (comma-separated):", 
                             placeholder = "Click on points in volcano plot to add gene names here", rows = 3),
               # textAreaInput(ns("typed_gene_list"), "Selected Genes (comma-separated):", 
               #               placeholder = "Type the gene names here", rows = 3),
               actionButton(ns("get_genes"), "plot genes"),
               actionButton(ns("clear_genes"), "clear")
        ),
        
        column(8, class = "plot-column",
               plotlyOutput(ns("gene_heatmap"), height = 'auto') # Plotly heatmap output area
               #textOutput(ns("selected_genes_text"))
        )
      )    
    )
  )
}


plot_server <- function(id, experiments, file_path, gene_data){#, molecule_type) { and experiments was user_info
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Create a unique source ID for this module instance
    plot_source <- paste0("volcano_", id)
    
    # Create isolated reactive values for this module instance
    local_selected <- reactiveValues(
      click = character(0),
      typed = character(0),
      gene_names = character(0),
      gene_ids = character(0),
      #is_active = TRUE,
      should_render_heatmap = FALSE,
      is_fresh_session = TRUE,  # Add flag for fresh session
      current_plot_ids = NULL,   # Add new reactive value for current plotted IDs
    )
    
    heatmap_data <- reactiveVal(NULL)
    
    #gene_data <- reactiveVal(gene_info)
    
    # Handle text input restoration
    # observeEvent(input$click_gene_list, {
    #   # Ignore if this is initial load or cleanup
    #   if (local_selected$is_fresh_session) {
    #     # Clear the text input if it has a value on fresh session
    #     if (!is.null(input$click_gene_list) && input$click_gene_list != "") {
    #       updateTextAreaInput(session, "click_gene_list", value = "")
    #     }
    #     return()
    #   }
    # }, ignoreInit = FALSE)  # We want to catch the initial value
    
    # Function to clean up the module state
    cleanup <- function() {
      print("=== Cleanup Event ===")
      
      # First disable rendering and active state
      local_selected$should_render_heatmap <- FALSE
      #local_selected$is_active <- FALSE
      local_selected$is_fresh_session <- TRUE  # Reset fresh session flag
      
      # Reset button tracking
      last_button_value(0)
      
      clear_genes()
      
      # Clear outputs
      output$gene_heatmap <- NULL
      output$volcano_plot <- NULL
      output$top_table <- NULL
      
      print("=== Cleanup Complete ===")
    }
    
    # # Reset fresh session flag on first real interaction
    # observeEvent(event_data("plotly_click", source = plot_source), {
    #   if (local_selected$is_fresh_session) {
    #     local_selected$is_fresh_session <- FALSE
    #   }
    # }, ignoreInit = TRUE)

    
    # Listen for cleanup message
    observeEvent(session$input[["plot_cleanup"]], {
      if (!is.null(session$input[["plot_cleanup"]]) && 
          session$input[["plot_cleanup"]]$id == id) {
        cleanup()
      }
    })
    
    # Initialize empty text areas for this instance
    # isolate({
    #   updateTextAreaInput(session, "click_gene_list", value = "")
    #   updateTextAreaInput(session, "typed_gene_list", value = "")
    # })

    # Add a reactive value to track initialization
    rv <- reactiveValues(initialized = FALSE)
    
    # Get all comparisons in the given experiment from the backend
    #all_comparison <- get_all_comparisons(user_info, exp_id)
    
    # Modify the comparison handling
    # observeEvent(all_comparison$out, {
    #   if (!is.null(all_comparison$out)) {
    #     updateSelectInput(session, "comparison",
    #                      choices = names(all_comparison$out), 
    #                      selected = names(all_comparison$out)[1])
    #     rv$initialized <- TRUE
    #   }
    # })
    # print(local_selected$is_active)
    # browser()
    # comparison_data <- reactive({
    #   browser()
    #   #req(local_selected$is_active) #, rv$initialized, input$comparison)
    #   #comp_id = all_comparison$out[input$comparison]
    #   #print(comp_id)
    #   browser()
    #   res = get_local_comparison_data(file_path)
    #   #molecule_table = query_output$molecule_table
    #   #print(colnames(res))
    #   browser()
    #   res$log10p = -log10(res$padj)
    #   res
    # })

    # sample condition information for the experiment
    # condition_data <- reactive({
    #   req(local_selected$is_active)
    #   res = get_sample_conditions(user_info, exp_id)
    #   print(dim(res$out))
    #   res$out
    # })
    # 
    #THIS DATA IS JUST THE DIFFERENT TYPE OF MOLECULE IDS (UNIPROT, NAME, ENSEMBLE)
    # molecule_data <- reactive({
    #   #req(local_selected$is_active)
    #   req(local_selected$genes)
    #   molecule_table <- get_heatmap_data(file_path, local_selected$genes)#exp_molecule(user_info, exp_id, molecule_type)$out
    #   #colnames(molecule_table)[2] = 'gene_name'
    #   #print(head(molecule_table))
    #   browser()
    #   molecule_table
    # })

    # plot_data = VOLCANO PLOT DATA with additional label column for up-regulated, down-regulated, or non-signficant gene
    plot_data <- reactive({
      res <- get_local_comparison_data(file_path)
      #res = merge(res, molecule_data(), by = 'id', all.x = TRUE)
      
      res$log10p <- log(res$padj, base=10) * -1

      # assign color dynamically, blue for down, red for up, gray for ns
      res$category = 'ns'
      res$category[res$log2FoldChange>=input$lfc_cutoff&res$log10p>= input$p_cutoff] = 'up'
      res$category[res$log2FoldChange<= -input$lfc_cutoff&res$log10p>= input$p_cutoff] = 'down'
      
      res$display_name <- ifelse(is.na(res$gene) | res$gene == "", res$gene_id, res$gene)
      #browser()
      res
    })
    
    # Render the interactive Plotly scatter plot with dynamic coloring
    color_dict = c("down" = "dodgerblue", "ns" = "gray", "up" = "darkorange")
    output$volcano_plot <- renderPlotly({
      withProgress(message = 'Retrieving data...', {
        req(plot_data())  # Ensure we have data before trying to plot
      })
      withProgress(message = 'Generating plot...', {
        p <- plot_ly(
          plot_data(), 
          x = ~log2FoldChange, y = ~log10p , key = ~paste(gene_id, gene, sep = "_"), text = ~display_name,# text = ~paste(gene_name, display_id, sep = ":"),
          type = 'scatter', mode = 'markers',
          marker = list(size = 10, opacity = 0.6),
          color = ~category,
          colors = color_dict,
          source = plot_source  # Use unique source ID
        ) %>%
        layout(
          title = list(
            text = paste("<b>Volcano Plot of ", input$comparison, "</b>"),
            font = list(size = 24, weight = "bold"), 
            y = 0.98  # Adjust title position
          ),
          margin = list(t = 50),  # Add top margin for title
          xaxis = list(
            title = list(
              text = "Log Fold Change",
              font = list(size = 18)  # X-axis title font size
            ),
            tickfont = list(size = 14)  # X-axis tick labels font size
          ),
          yaxis = list(
            title = list(
              text = "-log10(P)",
              font = list(size = 18)  # Y-axis title font size
            ),
            tickfont = list(size = 14)  # Y-axis tick labels font size
          ),
          legend = list(
            font = list(size = 16)  # Legend font size
          ),
          clickmode = 'event'
        )
        event_register(p,"plotly_click")
        
        observeEvent(event_data("plotly_click", source = plot_source), {
          #req(local_selected$is_active)
          click_data <- event_data("plotly_click", source = plot_source)
          if (!is.null(click_data)) {
            gene_id <- strsplit(click_data$key, "_")[[1]][1]
            gene_name <- strsplit(click_data$key, "_")[[1]][2]
            #browser()
            if(!gene_id %in% local_selected$gene_ids){
              local_selected$gene_names <- c(local_selected$gene_names, gene_name)
              local_selected$gene_names <- na.omit(local_selected$gene_names)
              local_selected$gene_ids <- c(local_selected$gene_ids, gene_id)
              local_selected$gene_ids <- na.omit(local_selected$gene_ids)
              print(paste("length of gene names", length(local_selected$gene_names)))
              print(paste("length of gene ids", length(local_selected$gene_ids)))
            }
            if(input$id_type == "gene_name")
              updateTextAreaInput(session, "genes", value = isolate(paste(local_selected$gene_names, collapse = ",")))
            else #use ensembl gene ID
              updateTextAreaInput(session, "genes", value = isolate(paste(local_selected$gene_ids, collapse = ",")))
          }
        })
        
        p
      })
      
    }) #END VOLCANO PLOT RENDER PLOTLY
    
    
    observeEvent(input$id_type, { #clear selected genes if new gene name type selected
      clear_genes()
    })
    
    #Use local_selected instead of selected in all other observers
    # observeEvent(input$click_gene_list, {
    #   modified_click_genes <- strsplit(input$click_gene_list, ",")[[1]]
    #   modified_click_genes <- trimws(modified_click_genes)
    #   modified_click_genes <- unique(modified_click_genes[modified_click_genes != ""])
    #   local_selected$click <- modified_click_genes
    # })
    # 
    # observeEvent(input$typed_gene_list, {
    #   typed_genes <- strsplit(input$typed_gene_list, ",")[[1]]
    #   typed_genes <- trimws(typed_genes)
    #   typed_genes <- unique(typed_genes[typed_genes != ""])
    #   local_selected$typed <- unique(typed_genes)
    # })

    # table_data <- reactive({
    #   sort_col_dict = c("log2FC", "log2FC", "pValue", "qValue")
    #   names(sort_col_dict) = c("Positive Log2FC", "Negative Log2FC", "P Value", "Adjusted P Value")
    #   sort_col = sort_col_dict[input$sort_by]
    #   #print(sort_col)
    #   cutoff_val = input$table_cutoff
    #   n_val = input$n_feature
    #   
    #   all_rows = plot_data()[,!colnames(plot_data()) %in% c('id','display_id')]
    #   #print(colnames(comparison_data))
    #   if (input$sort_by == "Negative Log2FC"){
    #     table_subset = all_rows[all_rows[,sort_col] <= -cutoff_val,]
    #     table_subset = table_subset[order(table_subset[,sort_col]),]
    #     
    #   }else if(input$sort_by == "Positive Log2FC"){
    #     table_subset = all_rows[all_rows[,sort_col] >= cutoff_val,]
    #     table_subset = table_subset[order(-table_subset[,sort_col]),]
    #   }else{
    #     table_subset = all_rows[all_rows[,sort_col] <= cutoff_val,]
    #     table_subset = table_subset[order(table_subset[,sort_col]),]
    #   }
    #   if(dim(table_subset)[1]>=n_val){
    #     table_subset = table_subset[0:n_val, ]
    #   }
    #   table_subset
    # })
    # 
    # output$top_table <- renderDT({
    #   datatable(table_data(),  
    #             extensions = 'Buttons',  # Add buttons extension
    #             options = list(
    #               pageLength = 10,
    #               dom = 'Bfrtip',  # Add 'B' for buttons
    #               buttons = list(  # Configure download buttons
    #                 list(
    #                   extend = 'csv',
    #                   filename = 'differential_analysis_results'
    #                 ),
    #                 list(
    #                   extend = 'excel',
    #                   filename = 'differential_analysis_results'
    #                 )
    #               ),
    #               autoWidth = TRUE,
    #               searchHighlight = TRUE,
    #               initComplete = JS(
    #                 "function(settings, json) {",
    #                 "$(this.api().table().container()).css({'font-size': '14px'});",
    #                 "}")
    #             ),
    #             rownames = F) %>%
    #     formatRound(columns = c("log2FC", "log10p", "log10q"), digits = 4) %>%
    #     formatRound(columns = c("pValue", "qValue"), digits = 6)
    # })

    # Query the molecule_data to get ids of both clicked and typed genes
    # selected_ids <- reactive({
    #   req(molecule_data()) #local_selected$is_active
    #   lookup = molecule_data()
    #   # print("Lookup table:")
    #   # print(head(lookup))
    #   # print(colnames(lookup))
    #   # print("ID type:")
    #   # print(input$id_type)
    # 
    #   # Handle clicked genes
    #   # clicked_ids = input$clicked_genes_list #lookup[as.character(lookup$display_id) %in% local_selected$click, ]
    #   # print("Clicked IDs:")
    #   # print(clicked_ids)
    # 
    #   # Handle typed genes - make case insensitive
    #   # typed_genes_lower = tolower(local_selected$typed)
    #   # lookup_col_lower = tolower(as.character(lookup[[input$id_type]]))
    #   # typed_ids = lookup[lookup_col_lower %in% typed_genes_lower, ]
    #   # print("Typed IDs:")
    #   # print(typed_ids)
    # 
    #   ids <- unique(c(as.character(clicked_ids$id), as.character(typed_ids$id)))
    #   selected_id_df <- lookup[as.character(lookup$id) %in% ids, c("id", "gene_name", "display_id")]
    #   selected_id_df$id <- as.character(selected_id_df$id)
    #   print("selected_id_df:")
    #   print(selected_id_df)
    # 
    #   list(ids = selected_id_df$id, gene_names = rownames(lookup))#selected_id_df$gene_name)#, display_id = selected_id_df$display_id)
    # })
    # Track previous button value
    #last_button_value <- reactiveVal(0)

    #Create the heatmap output only when get_genes is clicked (or when rendered?)
    observeEvent(input$get_genes, {
      #req(active_card()) #analysis results of file must already be generated
      req(local_selected$gene_ids) # there must be input in typed genes text input
      #browser()
      # print("=== Get Genes Event ===")
      # print(paste("Button value:", input$get_genes))
      #print(paste("Last button value:", last_button_value()))
      #print(paste("Should render heatmap:", local_selected$should_render_heatmap))

      # Skip if this is not a new button click
      # if (is.null(input$get_genes) || input$get_genes <= last_button_value()) {
      #   print("Skipping - not a new button click")
      #   return()
      # }

      #last_button_value(input$get_genes)

      #req(local_selected$is_active)
      #local_selected$genes = unique(c(local_selected$click, local_selected$typed))

      # Update current_plot_ids only when button is clicked
      # local_selected$current_plot_ids <- selected_ids()$ids

      # if (length(local_selected$current_plot_ids) == 0) {
      #   showNotification("No matching genes found", type = "warning")
      #   output$gene_heatmap <- NULL
      #   return()
      # }

      # Only set should_render_heatmap to TRUE on actual button clicks
      #local_selected$should_render_heatmap <- TRUE

      # Get sample value data
      withProgress(message = 'Retrieving molecule value data...', {
        if(grepl(file_path,"siNTC"))
          logcpm_file <- r"(experiments\siNTC_comparisons\logCPM_all_samples.txt)"
        else
          logcpm_file <- r"(experiments\ICP4_comparisons\logCPM_all_samples.txt)"
        df <- read.table(logcpm_file, header=TRUE)
        #browser()
        df_filtered <- df[rownames(df) %in% local_selected$gene_ids,]
        if(input$id_type == "gene_name"){
          new_rownames <- gene_data[rownames(df_filtered)]
          new_rownames <- make.unique(new_rownames, sep = "_") #make them unique, in case difference gene ids point to the same gene
          rownames(df_filtered) <- new_rownames
        }
        heatmap_data(df_filtered)
        #browser()
        
      })

      print("Creating heatmap...")
      output$gene_heatmap <- renderPlotly({
        withProgress(message = 'Preparing data...', {
          req(heatmap_data())
          #browser()
          # isolate({
          #   if(input$id_type == "gene_name"){ #adjust rownames of heatmap data based on user display preference
          #     heatmap_data({
          #       temp <- heatmap_data()
          #       # Get the new rownames and make them unique, in case difference gene ids point to the same gene
          #       new_rownames <- gene_data[rownames(temp)]
          #       new_rownames <- make.unique(new_rownames, sep = "_")
          #       browser()
          #       rownames(temp) <- new_rownames
          #       temp
          #     })
          #   }
          # })
          
          browser()
  
          # Ensure matrix structure when creating heatmap_data
          print(dim(heatmap_data()))
          print(colnames(heatmap_data()))
          print(rownames(heatmap_data()))
          # CHECKIF LINE BELOW IS NEEDED IF ONLY 1 SELECTED GENE
          # if (nrow(heatmap_data) == 1) {
          #   # Explicitly set dimensions for single-row case
          #   dim(heatmap_data) <- c(1, ncol(heatmap_data))
          #   rownames(heatmap_data) <- rownames(sample_data()$res)
          #   colnames(heatmap_data) <- colnames(sample_data()$res)
          # }
  
          data_range <- range(heatmap_data(), na.rm = TRUE, finite = TRUE)
          max_abs <- max(abs(data_range))
          if (min(data_range) < 0) {
            color_palette = colorRampPalette(c("dodgerblue2", "white", "darkred"))(256)
            breaks = seq(-max_abs, max_abs, length.out = 256)
          } else {
            color_palette = colorRampPalette(c("white", "indianred"))(256)
            breaks = seq(0, max_abs, length.out = 256)
          }

          data_type_title <- "logCPM" #sample_data()$data_type
          
          #browser()
          
  
          print("Final data range:")
          print(range(heatmap_data(), na.rm = TRUE, finite = TRUE))
          print("Breaks:")
          print(range(breaks))
  
          incProgress(1.0, detail = "Calculating dimensions...")
          num_rows <- nrow(heatmap_data())
          height_px <- max(400, num_rows * 40)
        }) #END OF GETTING HEATMAP DATA
        
        experiment_types <- sapply(strsplit(colnames(heatmap_data()), "_"), function(x) {
          if(length(x) > 1) {
            paste(x[-1], collapse = "_")
          } else {
            "Unknown"
          }
        })
        
        # Create color palette for experiment types
        unique_types <- unique(experiment_types)
        type_colors <- scales::hue_pal()(length(unique_types))
        names(type_colors) <- unique_types
        #col_colors <- type_colors[experiment_types]
        browser()
        
        col_annotation <- data.frame(
          Experiment = factor(experiment_types, levels = unique_types)
        )
        rownames(col_annotation) <- colnames(heatmap_data())

        withProgress(message = 'Generating heatmap...', {
          heatmaply(
            heatmap_data(),
            dendrogram = "both",
            col_side_colors = col_annotation,
            col_side_palette = type_colors,
            colors = color_palette,
            breaks = breaks,
            Rowv = T, Colv = T,
            grid_color = 'white',
            fontsize_row = 10,
            fontsize_col = 10,
            subplot_heights = c(0.1, 0.05, 0.85),  # Adjusted: dendrogram, annotation, heatmap#c(2/(num_rows + 2), num_rows/(2 + num_rows)),
            height = height_px,
            showticklabels = c(FALSE, TRUE),  # (columns, rows)
            limits = NULL,
            main = data_type_title,
            showlegend = TRUE,
            branches_lwd = 0.3,  # Control dendrogram line thickness
            label_names = c("Gene", "Sample", "Value") 
          ) %>%
          layout(
            xaxis = list(
              tickfont = list(size = 14, color = 'black'),
              showticklabels = FALSE  # Hide column labels
            ),
            yaxis = list(
              tickfont = list(size = 14, color = 'black'),
              showticklabels = FALSE   # Show row labels
            )
          )
        })
      }) #end of render_plotly
    }, ignoreInit = TRUE)  # end of observeEvent get_genes
    
    clear_genes <- function(){
      local_selected$gene_names <- character(0)
      local_selected$gene_ids <- character(0)
      updateTextAreaInput(session, "genes", value = "")
    }
    
      #Clear heatmap when clear_genes is clicked
    observeEvent(input$clear_genes, {
      #local_selected$should_render_heatmap <- FALSE
      clear_genes()
    })
    
    print(paste("Namespace:", session$ns("volcano_plot"))) 
    print('getting data')
    #browser()
    #print(paste('molecule_type: ', molecule_type))
    #Create reactiveVal to store the last selected IDs
    #last_selected_ids <- reactiveVal(NULL)
    
    # Update last_selected_ids when button is clicked
    # observeEvent(input$get_genes, {
    #   last_selected_ids(selected_ids())
    # })
    # # Create reactive that only depends on last_selected_ids
    # selected_ids_final <- reactive({
    #   last_selected_ids()
    # })
    
      # Return module interface
    return(list(
      cleanup = cleanup,
      local_selected = local_selected
    ))
  })
}