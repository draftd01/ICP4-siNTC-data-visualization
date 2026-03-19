library(bslib)

home_ui <- function(id) {
  ns <- NS(id)
  
  layout_column_wrap(
    width = 1,
    heights_equal = "row",
    
    card(
      #card_header("Summary"),
      
      card_body(
        style = "max-width: 900px; margin: 0 auto;",
        
        # Title
        h3(
          tags$a(href = "https://www.biorxiv.org/content/10.64898/2025.12.09.693233v1",
                 "Transcriptional consequences of herpes simplex virus 1 ICP4 inducible expression in uninfected cells"),
          style = "font-weight: 600;"
        ),
        
        # Authors
        p(
          class = "text-muted",
          style = "margin-bottom: 5px;",
          "Nora L Herzog¹,*,#, Sarah Keegan¹,#, Daniel Drafta¹, Rebecca Plessel², 
           Liam J Holt¹,³, Ian J Mohr², and Angus C Wilson²"
        ),
        
        # Affiliations
        p(
          style = "font-size: 0.9em; line-height: 1.4;",
          "¹ Institute for Systems Genetics, New York University Langone Health, New York, NY 10016, USA",
          tags$br(),
          "² Department of Microbiology, New York University School of Medicine, New York, NY 10016, USA",
          tags$br(),
          "³ Department of Biochemistry and Molecular Pharmacology, New York University Langone Health, New York, NY 10016, USA",
          tags$br(), tags$br(),
          "* Present address: Institute for Integrative Systems Biology, Universitat de València – CSIC, Spain",
          tags$br(),
          "# Nora L Herzog and Sarah Keegan contributed equally (alphabetical order)"
        ),
        
        # Correspondence
        p(
          style = "font-size: 0.9em;",
          strong("Correspondence: "),
          tags$a(href = "mailto:nora.herzog@uv.es",
                 "nora.herzog@uv.es")
        ),
        
        # Link to paper
        p(
          style = "font-size: 0.9em;",
          strong("Preprint: "),
          tags$a(href = "https://www.biorxiv.org/content/10.64898/2025.12.09.693233v1",
                 "https://www.biorxiv.org/content/10.64898/2025.12.09.693233v1")
        ),
        
        #hr(),
        
        # Abstract Header
        h5("Abstract", style = "margin-top: 20px;"),
        
        # Abstract text
        p(
          style = "text-align: justify; line-height: 1.6;",
          "Successful viral infections reflect the balanced outcome of a tightly regulated program of viral gene expression and 
            manipulation of the host cell to favor production of new infectious particles. The productive replication cycle of herpes 
            simplex virus 1 (HSV-1) is dependent on the immediate-early essential transcription factor ICP4, which positively or negatively
            regulates transcription of immediate-early, early, and late HSV-1 genes at different steps in the temporal cascade, largely 
            through sequence-specific binding to cis-acting regulatory elements. In contrast, direct regulation of host transcription 
            programming by ICP4 is less well understood. In this study we exogenously expressed doxycycline-inducible wild-type and 
            mutant ICP4 proteins in uninfected primary human fibroblasts and performed RNA-Seq to identify ICP4-driven changes to the 
            host transcriptome. Cross-referencing our findings to a published dataset of ICP4-dependent changes to the host transcriptome 
            in cells infected with HSV-1 provided validation for a subset of ICP4-regulated differentially-expressed genes. Furthermore, 
            mutations in ICP4 that disrupt interactions with the host transcriptional machinery or interfere with site-specific DNA 
            binding further modified this ICP4-driven remodeling of host transcription programming. Taken together, these data provide 
            a comprehensive transcriptomic analysis of ICP4-driven dysregulation of host gene expression in uninfected cells."
        )
      )
    )
  )
}

home_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    # Add server logic for the Home tab here
  })
}


