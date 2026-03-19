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
  #print(paste0("plot_ui id: ", id))
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
             #selectInput(ns("comparison"), "Select Comparison:", 
             #             choices = NULL),
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
      card_header(h5("All Sample Data")),
      fluidRow(
        column(4, 
               selectInput(ns("id_type"), "id Type:", choices = c(
                 "Gene Name" = "gene_name", #first is text in drop down, second is actual value submitted in form
                 "Ensembl Gene ID" = "ensembl_gene_id"
               )),
               selectizeInput(ns("genes"), "Selected Genes (type or select from plot)", choices=c("") , multiple = TRUE),
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


plot_server <- function(id, experiments, file_path, gene_id_to_name){#, molecule_type) { and experiments was user_info
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Create a unique source ID for this module instance
    plot_source <- paste0("volcano_", id)
    
    # Create isolated reactive values for this module instance
    local_selected <- reactiveValues(
      #click = character(0),
      #typed = character(0),
      gene_names = character(0),
      gene_ids = character(0),
      #is_active = TRUE,
      #should_render_heatmap = FALSE,
      #is_fresh_session = TRUE,  # Add flag for fresh session
      #current_plot_ids = NULL,   # Add new reactive value for current plotted IDs
    )
    clicked = reactiveVal(FALSE)
    
    heatmap_data <- reactiveVal(NULL)
    
    # Function to clean up the module state
    cleanup <- function() {
      print("=== Cleanup Event ===")
      
      # First disable rendering and active state
      #local_selected$should_render_heatmap <- FALSE
      #local_selected$is_active <- FALSE
      #local_selected$is_fresh_session <- TRUE  # Reset fresh session flag
      
      # Reset button tracking
      #last_button_value(0)
      
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
    
    observeEvent(input$id_type,{
      #print("setting up selectize")
      #browser()
      if (input$id_type == "gene_name")
        updateSelectizeInput(session,"genes",choices = unique(as.vector(gene_id_to_name)), server=TRUE)
      else
        updateSelectizeInput(session,"genes",choices = names(gene_id_to_name), server=TRUE)
    })

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

    # plot_data = VOLCANO PLOT DATA with additional label column for up-regulated, down-regulated, or non-signficant gene
    plot_data <- reactive({
      res <- get_local_comparison_data(file_path)
      
      res$log10p <- log(res$pvalue, base=10) * -1
      res$log10q <- log(res$padj, base=10) * -1

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
            text = paste("<b>Volcano Plot</b>"),# of ", input$comparison, "</b>"),
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
        
        p
      })
      
    }) #END VOLCANO PLOT RENDER PLOTLY
    
    observeEvent(event_data("plotly_click", source = plot_source), {
      clicked(TRUE)
      click_data <- event_data("plotly_click", source = plot_source)
      if (!is.null(click_data)) {
        gene_id <- strsplit(click_data$key, "_")[[1]][1]
        gene_name <- strsplit(click_data$key, "_")[[1]][2]
        #browser()
        current_genes <- input$genes
        #print("in click function")
        local_selected$gene_names <- union(local_selected$gene_names, gene_name)
        #print("got first set of gene names")
        local_selected$gene_names <- na.omit(local_selected$gene_names)
        local_selected$gene_ids <- union(local_selected$gene_ids, gene_id) #save gene ids for heatmap display later
        local_selected$gene_ids <- na.omit(local_selected$gene_ids)
        
        #browser()
        
        if(input$id_type == "gene_name"){
          #print("updating list with gene name")
          #print(union(current_genes,gene_name))
          updateSelectizeInput(session, "genes", choices = unique(as.vector(gene_id_to_name)), selected = union(current_genes,gene_name), server=T)
        }
        else {#use ensembl gene ID
          #print("updating list with gene id")
          #print(union(current_genes,gene_id))
          updateSelectizeInput(session, "genes", choices = names(gene_id_to_name), selected = union(current_genes,gene_id), server=T)
        }
        
      }
    })
    
    
    observeEvent(input$id_type, { #clear selected genes if new gene name type selected
      clear_genes()
    })

    table_data <- reactive({
      sort_col_dict = c("log2FoldChange", "log2FoldChange", "pvalue", "padj")
      names(sort_col_dict) = c("Positive Log2FC", "Negative Log2FC", "P Value", "Adjusted P Value")
      sort_col = sort_col_dict[input$sort_by]
      #print(sort_col)
      cutoff_val = input$table_cutoff
      n_val = input$n_feature
      
      #browser()
      all_rows = plot_data()[,!colnames(plot_data()) %in% c('SeqRun','baseMean','lfcSE','stat','category','display_name')]
      #print(colnames(comparison_data))
      if (input$sort_by == "Negative Log2FC"){
        table_subset = all_rows[all_rows[,sort_col] <= -cutoff_val,]
        table_subset = table_subset[order(table_subset[,sort_col]),]

      }else if(input$sort_by == "Positive Log2FC"){
        table_subset = all_rows[all_rows[,sort_col] >= cutoff_val,]
        table_subset = table_subset[order(-table_subset[,sort_col]),]
      }else{
        table_subset = all_rows[all_rows[,sort_col] <= cutoff_val,]
        table_subset = table_subset[order(table_subset[,sort_col]),]
      }
      if(dim(table_subset)[1]>=n_val){
        table_subset = table_subset[0:n_val, ]
      }
      table_subset
    })
    
    output$top_table <- renderDT({
      
      datatable(table_data(),
                extensions = 'Buttons',  # Add buttons extension
                options = list(
                  pageLength = 10,
                  dom = 'Bfrtip',  # Add 'B' for buttons
                  buttons = list(  # Configure download buttons
                    list(
                      extend = 'csv',
                      filename = 'differential_analysis_results'
                    ),
                    list(
                      extend = 'excel',
                      filename = 'differential_analysis_results'
                    )
                  ),
                  autoWidth = TRUE,
                  searchHighlight = TRUE,
                  initComplete = JS(
                    "function(settings, json) {",
                    "$(this.api().table().container()).css({'font-size': '14px'});",
                    "}")
                ),
                rownames = F) %>%
        formatRound(columns = c("log2FoldChange", "log10p", "log10q"), digits = 4) %>%
        formatRound(columns = c("pvalue", "padj"), digits = 6)
    })

    #Create the heatmap output only when get_genes is clicked (or when rendered?)
    observeEvent(input$get_genes, {
      req(input$genes)

      # Get sample value data
      withProgress(message = 'Retrieving molecule value data...', {
        if(grepl("mock",file_path))
          #logcpm_file <- r"(experiments\mock_comparisons\logCPM_all_samples.txt)"
          logcpm_file <- file.path("experiments","mock_comparisons","logCPM_all_samples.txt")
        else
          #logcpm_file <- r"(experiments\ICP4_comparisons\logCPM_all_samples.txt)"
          logcpm_file <- file.path("experiments","ICP4_comparisons","logCPM_all_samples.txt")
        df <- read.table(logcpm_file, header=TRUE) #columns are exp, rows are ensembl gene ids

        #rename colnames to match new filenames
        colnames(df) <- gsub("RR_LL", "mDBD", colnames(df))
        colnames(df) <- gsub("siNTC", "mock", colnames(df))
        
        experiment_cat <- sapply(strsplit(colnames(df), "_"), function(x) { #remove exp date from colnames, so date_exp -> exp
          paste(x[-1], collapse = "_")
        })
        exp_cat_mask <- sapply(experiment_cat, function(cat){ #T/F vector same length as exp_cat, not only true values
          grepl(cat, file_path)
        })
        
        #rename these columns after extracting experiment cateogories (which were based on old names)
        names(exp_cat_mask) <- gsub("n6", "dNTA", names(exp_cat_mask))
        names(exp_cat_mask) <- gsub("n208", "dCTA", names(exp_cat_mask))
        
        typed_gene_ids <- NULL
        if (input$id_type == "gene_name") {
          typed_genes <- setdiff(input$genes,local_selected$gene_names) #local_selected are genes that were clicked
          gene_ids <- gene_id_to_name[gene_id_to_name %in% typed_genes]
          typed_gene_ids <- names(gene_ids)
        }
        else{
          typed_gene_ids <- setdiff(input$genes,local_selected$gene_ids)
        }
        
        combined_gene_ids <- union(local_selected$gene_ids, typed_gene_ids)
        
        #browser()
        df_filtered <- df[rownames(df) %in% combined_gene_ids,exp_cat_mask]#local_selected$gene_ids, exp_cat_mask]
        
        if(input$id_type == "gene_name"){ #adjust rownames of heatmap data based on user display preference
          new_rownames <- gene_id_to_name[rownames(df_filtered)] #convert gene_ids to gene_names
          new_rownames <- make.unique(new_rownames, sep = "_") #make them unique, in case different gene ids point to the same gene
          rownames(df_filtered) <- new_rownames
        }
        new_colnames <- make.unique(names(exp_cat_mask[exp_cat_mask]),sep="-") #get only TRUE values from mask vector, and make names unique to be used as colnames
        colnames(df_filtered) <- new_colnames
        heatmap_data(df_filtered)
        
      })

      #print("Creating heatmap...")
      output$gene_heatmap <- renderPlotly({
        withProgress(message = 'Preparing data...', {
          req(heatmap_data())
  
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
  
          #print("Final data range:")
          #print(range(heatmap_data(), na.rm = TRUE, finite = TRUE))
          #print("Breaks:")
          #print(range(breaks))
  
          incProgress(1.0, detail = "Calculating dimensions...")
          num_rows <- nrow(heatmap_data())
          height_px <- max(400, num_rows * 40)
        }) #END OF GETTING HEATMAP DATA
        
        experiment_types <- sapply(strsplit(colnames(heatmap_data()), "-"), function(x) {
          x[1]
        })
        
        # Create color palette for experiment types
        unique_types <- unique(experiment_types)
        type_colors <- scales::hue_pal()(length(unique_types))

        names(type_colors) <- unique_types
        col_colors <- type_colors[experiment_types]
        
        col_annotation <- data.frame(
          Experiment = factor(experiment_types, levels = unique_types)
        )

        withProgress(message = 'Generating heatmap...', {
          heatmaply(
            heatmap_data(),
            Rowv = FALSE, 
            Colv = TRUE,
            dendrogram = "column",
            #dendrogram = "both",
            col_side_colors = col_annotation, #exp type groups
            col_side_palette = type_colors,  #color coding of exp types
            colors = color_palette,
            breaks = breaks,
            #Rowv = T, Colv = T,
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
            label_names = c("Gene", "Sample", "Value"),
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
      updateSelectizeInput(session,"genes",selected = c(""))
      #print("clearing genes")
    }
    
      #Clear heatmap when clear_genes is clicked
    observeEvent(input$clear_genes, {
      #local_selected$should_render_heatmap <- FALSE
      clear_genes()
    })
    
    observeEvent(input$genes, { #for when selected genes list gets altered (especially removed), update local_selected list of genes representing genes clicked from plot
      if (!clicked()){
        #print("local selected before")
        #print(local_selected$gene_names)
        #print(local_selected$gene_ids)
        #print("")
        if(input$id_type == "gene_name"){
          missing_gene_names <- setdiff(local_selected$gene_names,input$genes)
          missing_gene_ids <- names(gene_id_to_name[gene_id_to_name %in% missing_gene_names])
          local_selected$gene_names <- setdiff(local_selected$gene_names,missing_gene_names)
          local_selected$gene_ids <- setdiff(local_selected$gene_ids,missing_gene_ids)
        }else{
          missing_gene_ids <- setdiff(local_selected$gene_ids,input$genes)
          missing_gene_names <- as.vector(gene_id_to_name[missing_gene_ids])
          local_selected$gene_names <- setdiff(local_selected$gene_names,missing_gene_names)
          local_selected$gene_ids <- setdiff(local_selected$gene_ids,missing_gene_ids)
        }
        #print("local selected after")
        #print(local_selected$gene_names)
        #print(local_selected$gene_ids)
      }
      clicked(FALSE)
    })
    
      # Return module interface
    return(list(
      cleanup = cleanup,
      local_selected = local_selected
    ))
  })
}