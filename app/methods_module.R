library(bslib)

methods_ui <- function(id) {
  ns <- NS(id)
  layout_column_wrap(
    width = 1,
    heights_equal = "row",
    card(
      card_header ("Methods"),
      # -------------------------
      HTML("<b><u>Lentivirus Construction & Cell Treatment</u></b>"),
      p("Successful viral infections reflect the balanced outcome of a tightly regulated program of viral gene expression and 
        manipulation of the host cell environment to favor production of new infectious particles. The productive (lytic) replication
        cycle of herpes simplex virus 1 (HSV-1) is dependent on the essential transcription factor ICP4 encoded by one of five 
        immediate-early genes. At different steps in the HSV-1 temporal cascade, ICP4 either positively or negatively regulates
        transcription of immediate-early, early, and late HSV-1 genes, including its own, through sequence-specific binding to 
        cis-acting regulatory elements. Because of this, the direct regulatory consequences of ICP4 expression on the host 
        transcriptome are less well understood. In this study we used doxycycline-regulated lentiviral expression vectors to 
        inducibly express wild type and mutant ICP4 proteins in uninfected primary human fibroblasts and performed RNA-Seq to 
        identify ICP4-driven changes to the host transcriptome. Cross-referencing our findings to a published dataset of ICP4-dependent 
        changes to the host transcriptome in HSV-1 infected cells provided validation for a subset of differentially-expressed
        genes regulated by ICP4. Furthermore, disrupting the ICP4 DNA binding domain was sufficient to alter the cellular gene 
        transcriptional program responsive to ICP4. This indicates that the DNA-binding domain of ICP4, which is required for 
        site-specific DNA binding to the virus genome, may also regulate binding to the host genome. Together these data provide
        a comprehensive transcriptomic analysis of how wild type and mutant ICP4 expression impact cellular gene expression in uninfected cells."),
     # --------------------------
      HTML("<b><u>RNA-Seq Library Preparation & Analysis</u></b>"),
      HTML("<p>For RNA-Seq, extracted RNA samples were processed in the Genome Technology Center at NYU Langone Health. Sample quality 
        was determined using a Bioanalyzer before automated stranded RNA-Seq library prep with polyA selection on ~1 ug of RNA. 
        Sequencing was performed on NovaSeq X+ 10B 100 Cycle Flowcells. Paired-end sequencing reads (FASTQ format) of length 150 
        bp (n208/mDBD/ICP4 batch) or 50 bp (n6/ICP4 batch) were subjected to quality control using FastQC to assess read quality 
        metrics(39). Adapter/quality trimming and filtering were conducted with Trim Galore
        (<span><a href='https://github.com/FelixKrueger/TrimGalore' target='_blank'>https://github.com/FelixKrueger/TrimGalore</a></span>)(40)
        and minimum length and quality parameters adjusted after examining the FastQC report. Reads were aligned with STAR(41) to
        a modified GRCh38 reference genome that included the ICP4 WT or mutant sequence. STAR was run with default parameters with
        output in sorted BAM format (--outSAMtype BAM SortedByCoordinate) and gene-level quantification enabled (--quantMode GeneCounts).</p>"),

     p("Mock controls were sequenced as a separate batch (paired-end length 50 reads). QC/Trimming/Alignment was performed as 
       described above, aligning to unmodified GRCh38. Batch effect could not be modeled due to a lack of matching induced 
       genotypes from the sequencing run. However, examining PCA for each mutant with siNTC samples showed 86-89% variance on 
       PC1, the axis that separated siNTC from the mutant genotypes. ICP4 WT samples were sequenced with each mutant condition 
       and showed 70% variance accounted for on PC1, which separated siNTC samples from ICP4 WT samples. This was most likely 
       due to expression differences of WT ICP4 and/or technical differences from separate sequencing runs. While it is possible
       that these technical differences are detected in DE analysis, given that the largest variance was detected on the PC axis
       separating siNTC from ICP4-induced genotypes, we considered the differences found to be reflecting true biological signal."),
     
     p("Differential expression analysis was performed on all samples with a standard DESeq2 pipeline in R(42). Prior to DE 
       analysis, lowly expressed genes were removed by filtering to those where at least 10 reads were mapped in a minimum of 
       2 or 3 samples, depending on group size. We report raw p-values, q-values and log2FC for all remaining genes.  Linear 
       models included a term to account for batch effect, except when not possible (siNTC comparisons). GSEA pre-ranked(31) 
       was run for each DE comparison using the fgsea R package, ranking genes by sign(log2FC) x q-value."),
     
     HTML("<p>For comparisons between our datasets and existing ICP4 Chip-Seq in WT HSV-1 infection, ChipSeq data including Input 
       controls (FASTQ format) were downloaded from BioProject (<span><a href='https://www.ncbi.nlm.nih.gov/bioproject/PRJNA553563' target='_blank'>PRJNA553563</a></span>)
       and subjected to QC/Trimming/Alignment with FASTQC, Trim Galore and STAR (similar workflow as above), aligning to 
       GRCh38 with KOS (HSV-1) genome appended. Peaks were called with MACS2(43) for each of 2, 4 and 6 h time points, 
       producing narrowPeak files. The ChIPseeker R package was then used to annotate peak locations to genomic regions.
       Peaks were filtered to those found in Promotor and 5' UTR regions and then overlap with DEGs from RNASeq analysis by
       the nearest gene.</p>"),
     
     HTML("<p>For comparisons between our datasets and existing RNA-Seq in WT HSV-1 and n12 infection, Total RNA-Seq data including
       uninfected controls (FASTQ format) were downloaded from BioProject (<span><a href='https://www.ncbi.nlm.nih.gov/bioproject/PRJNA851702' target='_blank'>PRJNA851702</a></span>)
       and subjected to QC/Trimming/Alignment with FASTQC, Trim Galore and STAR (similar workflow as above), aligning to GRCh38
       with KOS (HSV-1) genome appended.</p>"),
     
     p("Over-representation analysis using the R package clusterProfiler was performed for pathway enrichment tests involving 
       unranked gene lists resulting from intersection of data sets(34).")
    ),
    card(
      card_header("References"),
      HTML("
        <p>8. DeLuca NA, Schaffer PA. 1988. Physical and functional domains of the herpes simplex virus transcriptional regulatory protein ICP4. J Virol 62:732–743.</p>
        <p>20. Shepard AA, Imbalzano AN, DeLuca NA. 1989. Separation of primary structural components conferring autoregulation, transactivation, and DNA-binding properties to the herpes simplex virus transcriptional regulatory protein ICP4. J Virol 63:3714–3728.</p>
        <p>26. Allen KE, Everett RD. 1997. Mutations which alter the DNA binding properties of the herpes simplex virus type 1 transactivating protein Vmw175 also affect its ability to support virus replication. J Gen Virol 78:2913–2922.</p>
        <p>39. Brown J, Pirrung M, McCue LA. 2017. FQC Dashboard: integrates FastQC results into a web-based, interactive, and extensible FASTQ quality control tool. Bioinformatics 33:3137–3139.</p>
        <p>40. Krueger F. TrimGalore: A wrapper around Cutadapt and FastQC to consistently apply adapter and quality trimming to FastQ files, with extra functionality for RRBS data. Github. <span><a href='https://github.com/FelixKrueger/TrimGalore' target='_blank'>https://github.com/FelixKrueger/TrimGalore</a></span>. Retrieved 18 August 2025.</p>
        <p>41. Dobin A, Davis CA, Schlesinger F, Drenkow J, Zaleski C, Jha S, Batut P, Chaisson M, Gingeras TR. 2013. STAR: ultrafast universal RNA-seq aligner. Bioinformatics 29:15–21.</p>
        <p>42. Love MI, Huber W, Anders S. 2014. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biol 15:550.</p>
        <p>31. Subramanian A, Tamayo P, Mootha VK, Mukherjee S, Ebert BL, Gillette MA, Paulovich A, Pomeroy SL, Golub TR, Lander ES, Mesirov JP. 2005. Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. Proc Natl Acad Sci U S A 102:15545–15550.</p>
        <p>43. Zhang Y, Liu T, Meyer CA, Eeckhoute J, Johnson DS, Bernstein BE, Nusbaum C, Myers RM, Brown M, Li W, Liu XS. 2008. Model-based analysis of ChIP-Seq (MACS). Genome Biol 9:R137.</p>
        <p>34. Yu G, Wang L-G, Han Y, He Q-Y. 2012. clusterProfiler: an R package for comparing biological themes among gene clusters. OMICS 16:284–287.</p>"
      )
    )
  )
}

methods_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    
  })
}