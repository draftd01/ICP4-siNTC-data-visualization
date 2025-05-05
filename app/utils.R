library(ggplot2)
library(plotly)
library(heatmaply)

# query functions
#url_base = "https://api.secrepedia.org/secrepediadb/"
url_base = "https://dev-api.secrepedia.org/secrepediadb/"

# Authentication function
with_auth <- function(operation_fn) {
  # Any request function always uses user_info
  function(user_info, ...){
   # make request with operation_fn and the access token
   make_request <- function(token){
    tryCatch({
      # the operation_fn should return a list with status, out, msg, handling known errors inside the operation_fn
      request_results <- operation_fn(token, ...)
      request_results
    },
    error = function(e){
      list(
        status = "error",  # Changed from status_code to status
        out = NULL,
        msg = sprintf("Unexpected error:%s",e$message)
      )
    })
   }

   # first try to make request with the current access token
   if (!is.null(user_info$access_token)){
    results <- make_request(user_info$access_token)
    
    # if unauthorized, try to refresh the token
    if (!is.null(results$status) && results$status == 401) {  # Added null check
      refresh_url <- paste0(url_base, "token/refresh/")
      refresh_response <- httr::POST(refresh_url, body = list(refresh = user_info$refresh_token), encode = "json")
      if (refresh_response$status_code == 200) {    # refresh success
          user_info$access_token <- httr::content(refresh_response)$access
          # try again with the new token
          results <- make_request(user_info$access_token)
      }
      else{
        # fail to refresh, urge user to login again
        results <- list(
          status = "error",
          out = NULL,
          msg = "Authentication failed, Please login again."
        )
      }
    }
   } 
   else {
    # no access token, send request without header
    results <- make_request(NULL)
   }
   return(results)
  }
}

show_exp_details <- with_auth(function(token, exp_id) {
  detail_url <- paste0(url_base, "experiment/", exp_id, "/")
  if (is.null(token)){
    detail_response <- httr::GET(detail_url, config = list())
  }
  else{
    detail_response <- httr::GET(detail_url, httr::add_headers(Authorization = sprintf("Bearer %s", token)))
  }
  if (detail_response$status_code == 200) {
    details <- httr::content(detail_response, "parsed", "application/json")
    status <- 200
    out <- details$description
    msg <- NULL
  }
  else {
    status <- detail_response$status_code
    out <- NULL
    msg <- "Error fetching experiment details."
  }
  return(list(status = status, out = out, msg = msg))
})

# Helper function to parse and get the molecule presence matrix
parse_molecule_presence <- function(molecule_presence_response){
  molecule_presence_data <- httr::content(molecule_presence_response, "parsed", "application/json")
  print("molecule_presence_data:")
  print(molecule_presence_data)
  molecule_ids  <- names(molecule_presence_data[[1]][[1]])
  exp_ids <- names(molecule_presence_data[[1]])
  print(exp_ids)

  if (length(molecule_ids) == 0){    # handle the case where no valid molecule ids are matched with the input
    res = NULL
  }
  else{
    print(exp_ids)
    print(molecule_ids)
    # Create list of dataframes with unique column names for presence
    res <- lapply(seq_along(molecule_presence_data[[1]]), function(i) {
      entry <- molecule_presence_data[[1]][[i]]
      df <- data.frame(id = names(entry),
                      presence = unlist(entry))
      # Use the exp_id as the column name
      colnames(df)[2] <- names(molecule_presence_data[[1]])[i]
      return(df)
    })
    
    # Merge all dataframes by 'id' column
    res <- Reduce(function(x, y) merge(x, y, by = "id", all = TRUE), res)
    
    # Set proper row and column names
    rownames(res) <- res$id
    res <- res[sort(molecule_ids), c('id',exp_ids)]
    print("res after sorting:")
    print(res)
    res$id <- NULL
    
    # For single experiment case, ensure matrix structure is correct
    if (length(exp_ids) == 1) {
      res <- matrix(res[[1]], nrow = 1)
      rownames(res) <- exp_ids
      colnames(res) <- sort(molecule_ids)

    } else {
      # Clean up row names before transpose
      res <- t(res)  # Transpose

    }
  }
  print(res)
  return(res)
}

# Helper function to query if a list of molecules are present in the list of experiments
search_molecules <- with_auth(function(token, id_type, molecule_ids, exp_ids){
  molecule_search_url <- paste0(url_base, "search-molecules/")
  payload <- list(id_type = id_type,
                  identifiers = as.list(molecule_ids),
                  experiment_ids = as.list(exp_ids))
  
  if (is.null(token)){
    molecule_search_response <- httr::POST(
      molecule_search_url, 
      body = payload,
      encode = "json",
      content_type("application/json")
    )
  } 
  else {
    molecule_search_response <- httr::POST(
      molecule_search_url, 
      body = payload,
      encode = "json",
      content_type("application/json"),
      httr::add_headers(Authorization = sprintf("Bearer %s", token))
    )
  }
  
  #print(molecule_search_response$status_code)
  if (molecule_search_response$status_code == 200){
    molecule_presence_data <- parse_molecule_presence(molecule_search_response)
    status <- 200
    out <- molecule_presence_data
    msg <- NULL
  }
  else if (molecule_search_response$status_code == 400){
    status <- 400
    error_content <- httr::content(molecule_search_response, "parsed")
    out <- NULL
    msg <- error_content$detail
  }
  else{
    status <- molecule_search_response$status_code
    out <- NULL
    msg <- "Error searching molecules."
  }
  return(list(status = status, out = out, msg = msg))
})

# Helper function to parse the comparison response to get the comparison label
comp_label <- function(comp){
  up = paste(sapply(comp$fc_up_conditions, function(cond){cond$short_name}), collapse = '-')
  down = paste(sapply(comp$fc_down_conditions, function(cond){cond$short_name}), collapse = '-')
  return(paste(up, down, sep = '_vs_'))
}

get_all_comparisons <- with_auth(function(token, exp_id){
  comp_url <- paste0(url_base, "comparison-data-sets/", exp_id, "/")
  if(is.null(token)){
    comp_response <- httr::GET(comp_url, config = list())
  }
  else{
    comp_response <- httr::GET(comp_url, httr::add_headers(Authorization = sprintf("Bearer %s", token)))
  }
  if (comp_response$status_code == 200){
    comps <- httr::content(comp_response, "parsed", "application/json")
    labels <- sapply(comps, comp_label)
    ids <- sapply(comps, function(comp){comp$id})
    out = ids
    names(out) = labels 
    msg = NULL
    status = 200
  }
  else{
    status <- comp_response$status_code
    out <- NULL
    msg <- "Error getting all comparisons from the given experiment."
  }
  return(list(status = status, out = out, msg = msg))
})


# Helper function to parse the comparison data response to get the comparison dataframe
parse_comparison_data <- function(comp_data_response, molecule_type){
  comp_data <- httr::content(comp_data_response, "parsed", "application/json")
  #print(comp_data)
  res <- data.frame(
    id = sapply(comp_data, function(entry){entry$molecule$id}),
    log2FC = sapply(comp_data, function(entry){if(is.null(entry$log2FC)) NA else entry$log2FC}),
    pValue = sapply(comp_data, function(entry){if(is.null(entry$pValue)) NA else entry$pValue}),
    qValue = sapply(comp_data, function(entry){if(is.null(entry$qValue)) NA else entry$qValue})
  )
  res$display_name = sapply(comp_data, function(entry){
    if(is.null(entry$molecule$display_name)) NA else entry$molecule$display_name})
  res$gene_biotype = sapply(comp_data, function(entry){
    if(is.null(entry$molecule$gene_biotype)) NA else entry$molecule$gene_biotype})
  res$gene_description = sapply(comp_data, function(entry){
    if(is.null(entry$molecule$gene_description)) NA else entry$molecule$gene_description})
  if (molecule_type == 'RNA'){
    # use ensembl_gene_id as display_id for RNA
    res$display_id = sapply(comp_data, function(entry){
      if(is.null(entry$molecule$ensembl_gene_id)) NA else entry$molecule$ensembl_gene_id})
  }
  else if (molecule_type == 'Protein'){
    res$display_id = sapply(comp_data, function(entry){
      if(is.null(entry$molecule$uniprot_id)) NA else entry$molecule$uniprot_id})
  }
  res = res[,c('id', 'display_id', 'display_name', 'gene_biotype', 'gene_description', 'log2FC', 'pValue', 'qValue')]

  return(res)
}

get_comparison_data <- with_auth(function(token, comp_id, molecule_type){
  comp_data_url <- paste0(url_base, "comparison-data/", comp_id, "/")
  if(is.null(token)){
    comp_data_response <- httr::GET(comp_data_url, config = list())
  }
  else{
    comp_data_response <- httr::GET(comp_data_url, httr::add_headers(Authorization = sprintf("Bearer %s", token)))
  }

  if (comp_data_response$status_code == 200){
    status <- 200
    out <- parse_comparison_data(comp_data_response, molecule_type)
    msg <- NULL
  }
  else if (comp_data_response$status_code == 404){
    status <- 404
    out <- NULL
    msg <- "No ComparisonData matches the given query."
  }
  else{
    status <- comp_data_response$status_code
    out <- NULL
    msg <- "Error getting comparison data."
  }
  return(list(status = status, out = out, msg = msg))
})

# Helper function to parse the sample conditions from response
parse_sample_conditions <- function(sample_response, molecule_type){
  sample_data <- httr::content(sample_response, "parsed", "application/json")
  sample_ids = sapply(sample_data, function(entry){entry$sample$id})
  sample_names = sapply(sample_data, function(entry){entry$sample$name})
  replicate_ids = sapply(sample_data, function(entry){entry$sample$replicate_number})
  library_size = sapply(sample_data, function(entry){if(is.null(entry$library_size)) NA else entry$library_size})
  
  n_conditions = sapply(sample_data, function(entry){length(entry$sample$condition)})
  n_conditions = max(unique(n_conditions))
  condition_categories = rep(NA, n_conditions)
  condition_names = rep(NA, n_conditions)
  condition_values = list()
  for (i in 1:n_conditions){
    condition_categories[i] = unique(sapply(sample_data, function(entry){entry$sample$condition[[i]]$category_name}))
    condition_values[[i]] = sapply(sample_data, function(entry){entry$sample$condition[[i]]$name})
    #print(condition_values[[i]])
  }
  condition_df = data.frame(condition_values)
  colnames(condition_df) = condition_categories
  print(paste("condition categories:", paste(condition_categories, collapse=",")))
  condition_df$sample_name = sample_names
  condition_df$replicate_id = replicate_ids
  condition_df$library_size = library_size    # library size for RNA data, used for logcpm calculation when data unit isread_counts
  print(paste("sample ids:", paste(sample_ids, collapse=",")))
  rownames(condition_df) = sample_ids
  return(condition_df)
}

# Helper function to get sample conditions for a given experiment id
get_sample_conditions <- with_auth(function(token, exp_id){
  sample_url <- paste0(url_base, "sample-data-sets/", exp_id, "/")
  if(is.null(token)){
    sample_response <- httr::GET(sample_url, config = list())
  }
  else{
    sample_response <- httr::GET(sample_url, httr::add_headers(Authorization = sprintf("Bearer %s", token)))
  }
  
  if (sample_response$status_code == 200){
    status <- 200
    out <- parse_sample_conditions(sample_response)
    msg <- NULL
  }
  else{
    status <- sample_response$status_code
    out <- NULL
    msg <- "Error getting sample conditions."
  }
  return(list(status = status, out = out, msg = msg))
})


# Helper function to parse the molecule information
parse_exp_molecule <- function(exp_molecule_response, molecule_type){
    exp_molecule_data <- httr::content(exp_molecule_response, "parsed", "application/json")
    res = data.frame(id = sapply(exp_molecule_data, function(entry){entry$id}),
                     display_name = sapply(exp_molecule_data, function(entry){entry$display_name}))
    if (molecule_type == 'RNA'){
      res$gene_name = sapply(exp_molecule_data, function(entry){if (is.null(entry$name)) NA else entry$name})
      res$ensembl_gene_id = sapply(exp_molecule_data, function(entry){if (is.null(entry$ensembl_id)) NA else entry$ensembl_id})
      res$display_id = res$ensembl_gene_id    # this is for showing unique ids in the input field of the intensity heatmap
      
    }
    else if (molecule_type == 'Protein'){
      res$gene_name = sapply(exp_molecule_data, function(entry){
        if (is.null(entry$linked_gene$name)) NA else entry$linked_gene$name})
      res$ensemble_peptide_id = sapply(exp_molecule_data, function(entry){
        if (is.null(entry$ensembl_id)) NA else entry$ensembl_id})
      res$uniprot_id = sapply(exp_molecule_data, function(entry){
        if (is.null(entry$uniprot_id)) NA else entry$uniprot_id})
      res$ensembl_gene_id = sapply(exp_molecule_data, function(entry){
        if (is.null(entry$linked_gene$ensembl_id)) NA else entry$linked_gene$ensembl_id})
      res$display_id = res$uniprot_id    # this is for showing unique ids in the input field of the intensity heatmap
    }
    return(res)
}

# Helper function to get all molecule information for a given experiment id
get_exp_molecule <- with_auth(function(token, exp_id, molecule_type){
  url <- paste0(url_base, "experiment-molecules/", exp_id, "/" )
  if(is.null(token)){
    exp_molecule_response <- httr::GET(url, config = list())
  }
  else{
    exp_molecule_response <- httr::GET(url, httr::add_headers(Authorization = sprintf("Bearer %s", token)))
  }
  if (exp_molecule_response$status_code == 200){
    status <- 200
    out <- parse_exp_molecule(exp_molecule_response, molecule_type)
    msg <- NULL
  }
  else if (exp_molecule_response$status_code == 404){
    status <- 404
    out <- NULL
    msg <- "No Experiment matches the given query."
  }
  else{
    status <- exp_molecule_response$status_code
    out <- NULL
    msg <- "Error getting molecule information."
  }
  return(list(status = status, out = out, msg = msg))
})


# Helper function to parse the sample-molecule value
parse_sample_molecule_values <- function(sample_data_response, molecule_type){
  sample_value_data <- httr::content(sample_data_response, "parsed", "application/json")
  #print(sample_value_data[[1]])
  data_type = NA
  if (molecule_type == 'RNA'){
    if (!is.null(sample_value_data[[1]]$logcpm)){
      print('logcpm')
      res = data.frame(id = sapply(sample_value_data, function(entry){entry$molecule$id}),
                       dataset_id = sapply(sample_value_data, function(entry){
                        if(is.null(entry$data_set_id)) NA else entry$data_set_id}),
                       logcpm = sapply(sample_value_data, function(entry){
                        if(is.null(entry$logcpm)) NA else entry$logcpm}))
      data_type = 'logcpm'
    }
    else{
      print('readcounts')
      res = data.frame(id = sapply(sample_value_data, function(entry){entry$molecule$id}),
                       dataset_id = sapply(sample_value_data, function(entry){
                        if(is.null(entry$data_set_id)) NA else entry$data_set_id}),
                       read_counts = sapply(sample_value_data, function(entry){
                        if(is.null(entry$read_counts)) NA else entry$read_counts}))
      data_type = 'read_counts'
    }
  }
  else if (molecule_type == 'Protein'){
    if (!is.null(sample_value_data[[1]]$raw_intensity)){
      print('raw intensity')
      res = data.frame(id = sapply(sample_value_data, function(entry){entry$molecule$id}),
                       dataset_id = sapply(sample_value_data, function(entry){
                        if(is.null(entry$data_set_id)) NA else entry$data_set_id}),
                       raw_intensity = sapply(sample_value_data, function(entry){
                        if(is.null(entry$raw_intensity)) NA else entry$raw_intensity}))
      data_type = 'raw_intensity'
    }
    else{
      print('log intensity')
      res = data.frame(id = sapply(sample_value_data, function(entry){entry$molecule$id}),
                       dataset_id = sapply(sample_value_data, function(entry){
                        if(is.null(entry$data_set_id)) NA else entry$data_set_id}),
                       log_intensity = sapply(sample_value_data, function(entry){
                        if(is.null(entry$log_intensity)) NA else entry$log_intensity}))
      data_type = 'log_intensity'
    }
  }
  
  res_wide = reshape(res, 
                   direction = "wide",
                   idvar = "id",
                   timevar = "dataset_id",
                   v.names = data_type)
  print(paste("res_wide dimensions:", paste(dim(res_wide), collapse="x")))
  colnames(res_wide) = gsub(paste0(data_type, "."), "", colnames(res_wide))
  return(list(res = res_wide, data_type = data_type))
}

# Helper function to get molecule values for a list of SampleDataSet ids and a list of molecule ids
get_exp_molecule_values <- with_auth(function(token, exp_id, molecule_type, molecule_ids){
  exp_molecule_data_url <- paste0(url_base, "experiment-sample-data/", exp_id, "/")
  # Ensure molecule_ids is always a list with at least one element
  molecule_ids <- as.numeric(molecule_ids)
  payload <- list(
    identifiers = I(molecule_ids)  # Use I() to prevent single element from being unboxed
  )

  if(is.null(token)){
    exp_molecule_data_response <- httr::POST(
      exp_molecule_data_url, 
      body = payload,
      encode = "json",
      content_type("application/json")
    )
  }
  else{
    exp_molecule_data_response <- httr::POST(
      exp_molecule_data_url, 
      body = payload,
      encode = "json",
      content_type("application/json"),
      httr::add_headers(Authorization = sprintf("Bearer %s", token))
    )
  }
  
  if (exp_molecule_data_response$status_code == 200) {
    status <- 200
    out <- parse_sample_molecule_values(exp_molecule_data_response, molecule_type)
    msg <- NULL
  }

  else{
    status <- exp_molecule_data_response$status_code
    out <- NULL
    msg <- rawToChar(exp_molecule_data_response$content)
  }
  return(list(status = status, out = out, msg = msg))
})


# Helper function to get data for secretion prediction
get_prediction_data <- function(id_type, ids){
  # send request to backend to get prediction data
  payload <- list(
    id_type = id_type,
    identifiers = as.list(ids)  # Convert vector to list to ensure proper JSON array serialization
  )
  prediction_url <- paste0(url_base, "secretion-prediction-data/")
  prediction_response <- httr::POST(prediction_url, body = payload, encode = "json")
  #print(prediction_response$status_code)
  query_data <- httr::content(prediction_response, "parsed", "application/json")
  return(query_data)
}

# Helper function to parse the prediction data response to get the prediction dataframe
parse_prediction_data <- function(query_data){
  # Add error checking for empty query data
  if(is.null(query_data) || length(query_data) == 0) {
    return(NULL)
  }
   
  # Safely create gene_res dataframe with error checking
  gene_res <- tryCatch({
    data.frame(
      name = sapply(query_data, function(entry){
        if(is.null(entry) || is.null(entry$display_name)) NA else entry$display_name
      }),
      ensembl_peptide_id = sapply(query_data, function(entry){
        if(is.null(entry) || is.null(entry$ensembl_id)) NA else entry$ensembl_id
      }),
      uniprot_id = sapply(query_data, function(entry){
        if(is.null(entry) || is.null(entry$uniprot_id)) NA else entry$uniprot_id
      }),
      ensembl_gene_id = sapply(query_data, function(entry){
        if(is.null(entry) || is.null(entry$linked_gene) || is.null(entry$linked_gene$ensembl_id)) NA else entry$linked_gene$ensembl_id
      }),
      organism = sapply(query_data, function(entry){
        if(is.null(entry) || is.null(entry$prediction_data) || length(entry$prediction_data) == 0) NA 
        else if(is.null(entry$prediction_data[[1]]$organism)) NA 
        else entry$prediction_data[[1]]$organism
      }),
      description = sapply(query_data, function(entry){
        if(is.null(entry) || is.null(entry$linked_gene) || is.null(entry$linked_gene$description)) NA else entry$linked_gene$description
      }),
      bio_type = sapply(query_data, function(entry){
        if(is.null(entry) || is.null(entry$linked_gene) || is.null(entry$linked_gene$biotype)) NA else entry$linked_gene$biotype
      })
    )
  }, error = function(e) {
    print(paste("Error creating gene_res:", e$message))
    return(NULL)
  })

  if(is.null(gene_res)) return(NULL)

  parse_algorithm_data <- function(prediction_data){
    if(is.null(prediction_data) || length(prediction_data) == 0) {
      text_algorithms = c('TOPCONS','Phobius', 'Predisi', 'SignalP','Outcyte', 'DeepSig')
      value_algorithms = 'Score'
      return(list(
        text_out = setNames(rep(NA, length(text_algorithms)), text_algorithms),
        value_out = setNames(rep(NA, length(value_algorithms)), value_algorithms)
      ))
    }

    tryCatch({
      algorithm_data <- data.frame(
        algorithm = sapply(prediction_data, function(entry){
          if(is.null(entry) || is.null(entry$algorithm)) NA else entry$algorithm
        }),
        value = sapply(prediction_data, function(entry){
          if(is.null(entry) || is.null(entry$value)) NA else entry$value
        }),
        text_value = sapply(prediction_data, function(entry){
          if(is.null(entry) || is.null(entry$text_value)) NA else entry$text_value
        })
      )

      text_algorithms = c('TOPCONS','Phobius', 'Predisi', 'SignalP','Outcyte', 'DeepSig')
      value_algorithms = 'Score'
      
      # Handle missing algorithms
      all_algorithms <- unique(c(text_algorithms, value_algorithms))
      missing_algorithms <- setdiff(all_algorithms, algorithm_data$algorithm)
      
      if(length(missing_algorithms) > 0) {
        missing_df <- data.frame(
          algorithm = missing_algorithms,
          value = NA,
          text_value = NA
        )
        algorithm_data <- rbind(algorithm_data, missing_df)
      }
      
      rownames(algorithm_data) <- algorithm_data$algorithm
      
      text_out <- setNames(rep(NA, length(text_algorithms)), text_algorithms)
      text_out[intersect(names(text_out), algorithm_data$algorithm)] <- 
        algorithm_data[intersect(names(text_out), algorithm_data$algorithm), 'text_value']
      
      value_out <- setNames(rep(NA, length(value_algorithms)), value_algorithms)
      value_out[intersect(names(value_out), algorithm_data$algorithm)] <- 
        algorithm_data[intersect(names(value_out), algorithm_data$algorithm), 'value']
      
      return(list(text_out = text_out, value_out = value_out))
    }, error = function(e) {
      print(paste("Error in parse_algorithm_data:", e$message))
      text_algorithms = c('TOPCONS','Phobius', 'Predisi', 'SignalP','Outcyte', 'DeepSig')
      value_algorithms = 'Score'
      return(list(
        text_out = setNames(rep(NA, length(text_algorithms)), text_algorithms),
        value_out = setNames(rep(NA, length(value_algorithms)), value_algorithms)
      ))
    })
  }
  
  # Process algorithm results with error handling
  algorithm_res <- lapply(query_data, function(entry) {
    tryCatch({
      if(is.null(entry) || is.null(entry$prediction_data)) {
        return(parse_algorithm_data(NULL))
      }
      parse_algorithm_data(entry$prediction_data)
    }, error = function(e) {
      print(paste("Error processing entry:", e$message))
      return(parse_algorithm_data(NULL))
    })
  })

  # Safely create dataframes
  tryCatch({
    text_df <- do.call(rbind, lapply(algorithm_res, function(x) {
      if(is.null(x) || is.null(x$text_out)) {
        return(data.frame(t(setNames(rep(NA, 6), c('TOPCONS','Phobius', 'Predisi', 'SignalP','Outcyte', 'DeepSig')))))
      }
      data.frame(t(x$text_out))
    }))
    
    value_df <- do.call(rbind, lapply(algorithm_res, function(x) {
      if(is.null(x) || is.null(x$value_out)) {
        return(data.frame(Score = NA))
      }
      data.frame(t(x$value_out))
    }))
    
    # Combine all dataframes
    final_df <- cbind(gene_res, text_df, value_df)
    return(final_df)
    
  }, error = function(e) {
    print(paste("Error creating final dataframe:", e$message))
    return(NULL)
  })
}
