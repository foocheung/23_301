#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(shiny)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(ggplot2)
  library(ggrepel)
  library(DT)
  library(forcats)
  library(DiagrammeR)
})

# =============================================================================
# CONFIG
# =============================================================================

#pbmc_out_dir <- "03_GSEA"
pbmc_out_dir <- "deg_gsea_pooled_vs_within_pool"

disease_folder   <- "all_apeced_vs_hc"
disease_contrast <- "MildUntreated_vs_HC_subject_level"

treatment_folder   <- "paired_cross_pool"
treatment_contrast <- "MildTreated_vs_MildUntreated_paired"

databases <- c("GO_BP", "Reactome", "Hallmark")

# =============================================================================
# EMBEDDED ANALYSIS DATA FLOW
# =============================================================================

analysis_flow_gv <- paste(c(
  "digraph Dynamic_Graph {",
  "  graph [layout = dot, rankdir = TB, nodesep = 0.8, ranksep = 1.0, splines = ortho, fontsize = 36, label = \"PBMC Treatment-Associated Pathway Filtering\", labelloc = t]",
  "  node [fontname = \"Arial\", fontsize = 34]",
  "  edge [fontname = \"Arial\", fontsize = 32]",
  "  gsea [label = \"GSEA results",
  "Treated vs Untreated + Untreated vs Healthy\", shape = box, style = filled, fillcolor = lightcoral]",
  "  app [label = \"Custom PBMC pathway-filtering app\", shape = box, style = filled, fillcolor = lightyellow]",
  "  treat [label = \"Treatment response",
  "Treatment FDR < 0.05\", shape = box, style = filled, fillcolor = pink]",
  "  recovery [label = \"Remove recovery",
  "Opposite Disease/Treatment NES directions\", shape = box, style = filled, fillcolor = pink]",
  "  disease [label = \"Select non-recovery candidates",
  "Disease FDR >= 0.05\", shape = box, style = filled, fillcolor = pink]",
  "  q1 [label = \"QUESTION 1 ANSWERED",
  "Pathways changed with treatment,",
  "not significant in untreated vs healthy\", shape = note, style = filled, fillcolor = lightyellow]",
  "  priority [label = \"Prioritize strongest effects",
  "|Treatment NES| >= 1.5",
  "|Disease NES| <= 1.0\", shape = box, style = filled, fillcolor = pink]",
  "  crossdb [label = \"Assess biological-program support",
  "across pathway databases\", shape = box, style = filled, fillcolor = pink]",
  "  q2 [label = \"QUESTION 2 ANSWERED",
  "Strongest pathways/programs",
  "meeting the criteria\", shape = note, style = filled, fillcolor = lightyellow]",
  "  interpret [label = \"Biological interpretation",
  "Identify non-recovery treatment effects",
  "that may reflect broader JAK/STAT inhibition\", shape = box, style = filled, fillcolor = lightblue]",
  "  q3 [label = \"QUESTION 3 ANSWERED",
  "Candidate broader pharmacologic effects",
  "that may be avoidable with more targeted therapy\", shape = note, style = filled, fillcolor = lightyellow]",
  "  gsea -> app []",
  "  app -> treat []",
  "  treat -> recovery []",
  "  recovery -> disease []",
  "  disease -> q1 []",
  "  disease -> priority []",
  "  priority -> crossdb []",
  "  crossdb -> q2 []",
  "  crossdb -> interpret []",
  "  interpret -> q3 []",
  "}"
), collapse = "\n")

# =============================================================================
# HELPERS
# =============================================================================

num_col <- function(x) suppressWarnings(as.numeric(x))

# Descriptive NES extremeness percentile.
# Calculated from |NES| against ALL tested pathways within the same
# cell type + database + comparison BEFORE Step 1/2/3 filtering.
# 99 means the pathway has a larger |NES| than ~99% of pathways in that stratum.
nes_extremeness_percentile <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  ok <- is.finite(x)
  out <- rep(NA_real_, length(x))
  n_ok <- sum(ok)

  if (n_ok == 0) return(out)
  if (n_ok == 1) {
    out[ok] <- 100
    return(out)
  }

  out[ok] <- 100 * dplyr::percent_rank(abs(x[ok]))
  out
}

get_pbmc_celltypes <- function(folder) {
  d <- file.path(pbmc_out_dir, folder, "gsea_by_celltype")
  if (!dir.exists(d)) stop("Missing folder: ", d)
  basename(list.dirs(d, recursive = FALSE, full.names = TRUE))
}

find_pbmc_gsea_file <- function(folder, celltype_safe, contrast, database) {
  file.path(
    pbmc_out_dir, folder, "gsea_by_celltype", celltype_safe,
    paste0(celltype_safe, "_", contrast, "_", database, "_results.csv")
  )
}

celltype_from_safe <- function(x) str_replace_all(x, "_", " ")

clean_label <- function(x) {
  x %>%
    str_remove("^GOBP_") %>%
    str_remove("^REACTOME_") %>%
    str_remove("^HALLMARK_") %>%
    str_replace_all("_", " ") %>%
    str_squish() %>%
    str_to_sentence()
}

assign_category <- function(x) {
  z <- str_to_lower(x)

  case_when(
    # Type II interferon / IFN-gamma
    str_detect(z, paste0(
      "interferon gamma|ifn[-_ ]?gamma|ifng|",
      "type ii interferon|type 2 interferon|",
      "response to interferon gamma|",
      "interferon-gamma"
    )) ~ "IFN-gamma / type II interferon signaling",

    # Type I interferon / IFN-alpha-beta
    str_detect(z, paste0(
      "interferon alpha|ifn[-_ ]?alpha|ifna|",
      "interferon beta|ifn[-_ ]?beta|ifnb|",
      "type i interferon|type 1 interferon|",
      "response to type i interferon|",
      "response to interferon alpha|",
      "response to interferon beta"
    )) ~ "IFN-alpha/beta / type I interferon signaling",

    # Broad/shared interferon parent terms that do not specify type I or II.
    # Example: Reactome "Interferon Signaling".
    str_detect(z, paste0(
      "(^|[-_ ])interferon[-_ ]+signaling($|[-_ ])|",
      "(^|[-_ ])interferon[-_ ]+signalling($|[-_ ])"
    )) ~ "Broad / shared interferon signaling",

    # Antiviral programs that are not explicitly type I or type II IFN
    str_detect(z, paste0(
      "response to virus|defense response to virus|",
      "antiviral|viral process|viral genome|",
      "virus infection"
    )) ~ "Antiviral response",

    # Other cytokine/JAK-STAT signaling
    str_detect(z, paste0(
      "jak|stat|cytokine|interleukin|il[-_ ]?[0-9]+|",
      "growth factor signaling|gp130|",
      "signal transduction by cytokine"
    )) ~ "Other cytokine / JAK-STAT signaling",

    str_detect(z, paste0(
      "translation|ribosom|40s|60s|peptide chain|",
      "nonsense mediated decay|nmd|no-go|non-stop|",
      "rqt|eif2|selenocysteine|selenoamino|",
      "protein targeting to membrane|srp-dependent"
    )) ~ "Translation / ribosome / RNA quality control",

    str_detect(z, paste0(
      "cell cycle|mitotic|mitosis|m phase|s phase|",
      "g1/s|g2/m|chromosome segregation|chromatid|",
      "spindle checkpoint|spindle assembly|dna replication|",
      "dna synthesis|e2f|nucleosome assembly|centromere|cenpa"
    )) ~ "Cell cycle / proliferation",

    str_detect(z, paste0(
      "antigen processing|antigen presentation|",
      "mhc|major histocompatibility|peptide antigen"
    )) ~ "Antigen processing / presentation",

    str_detect(z, paste0(
      "immune response|immune system|immunoregulatory|",
      "lymphocyte|t cell|b cell|monocyte|macrophage|",
      "dendritic|neutrophil|chemotaxis|complement|",
      "response to bacterium|defense response to bacterium"
    )) ~ "Other immune / innate signaling",

    str_detect(z, paste0(
      "oxidative phosphorylation|electron transport|",
      "cellular respiration|mitochond|glycolysis|",
      "amino acid|metabolism|metabolic|",
      "unfolded protein|endoplasmic reticulum stress|",
      "chaperone|protein folding"
    )) ~ "Metabolism / cellular stress",

    TRUE ~ "Other"
  )
}

read_one <- function(celltype_safe, celltype, database, folder, contrast, prefix) {
  f <- find_pbmc_gsea_file(folder, celltype_safe, contrast, database)
  if (!file.exists(f)) return(NULL)

  x <- suppressWarnings(fread(f, fill = TRUE))
  needed <- c("ID", "Description", "NES", "p.adjust")
  if (nrow(x) == 0 || !all(needed %in% names(x)) || "reason" %in% names(x)) {
    return(NULL)
  }

  x %>%
    transmute(
      celltype,
      celltype_safe,
      database,
      feature_id    = as.character(ID),
      feature_label = as.character(Description),
      NES           = num_col(NES),
      FDR           = num_col(p.adjust)
    ) %>%
    rename_with(~ paste0(prefix, "_", .x), c("NES", "FDR"))
}

# =============================================================================
# LOAD ONCE AT APP START
# =============================================================================

celltype_key <- tibble(
  celltype_safe = get_pbmc_celltypes(treatment_folder)
) %>%
  mutate(celltype = celltype_from_safe(celltype_safe))

treatment_all <- tidyr::expand_grid(celltype_key, database = databases) %>%
  pmap_dfr(function(celltype_safe, celltype, database) {
    read_one(
      celltype_safe, celltype, database,
      treatment_folder, treatment_contrast, "treatment"
    )
  })

disease_all <- tidyr::expand_grid(celltype_key, database = databases) %>%
  pmap_dfr(function(celltype_safe, celltype, database) {
    read_one(
      celltype_safe, celltype, database,
      disease_folder, disease_contrast, "disease"
    )
  })

# -------------------------------------------------------------------------
# NES EXTREMENESS CONTEXT
# -------------------------------------------------------------------------
# IMPORTANT: these percentiles are computed BEFORE any pathway filtering.
# Each pathway is compared only with all other tested pathways from the SAME
# cell type, pathway database, and comparison. This avoids comparing raw NES
# values across databases with different NES distributions.
treatment_all <- treatment_all %>%
  group_by(celltype, database) %>%
  mutate(
    treatment_abs_NES_percentile =
      nes_extremeness_percentile(treatment_NES)
  ) %>%
  ungroup()

disease_all <- disease_all %>%
  group_by(celltype, database) %>%
  mutate(
    disease_abs_NES_percentile =
      nes_extremeness_percentile(disease_NES)
  ) %>%
  ungroup()

# Fail early if duplicate pathway keys would multiply the join.
dup_t <- treatment_all %>%
  count(celltype, celltype_safe, database, feature_id) %>%
  filter(n > 1)

dup_d <- disease_all %>%
  count(celltype, celltype_safe, database, feature_id) %>%
  filter(n > 1)

if (nrow(dup_t) > 0 || nrow(dup_d) > 0) {
  stop("Duplicate pathway keys detected in treatment or disease inputs.")
}

joined <- treatment_all %>%
  left_join(
    disease_all %>%
      select(
        celltype, celltype_safe, database, feature_id,
        disease_NES, disease_FDR, disease_abs_NES_percentile
      ),
    by = c("celltype", "celltype_safe", "database", "feature_id")
  ) %>%
  mutate(
    abs_treatment_NES = abs(treatment_NES),
    abs_disease_NES   = abs(disease_NES),
    extremeness_gap_percentile =
      treatment_abs_NES_percentile - disease_abs_NES_percentile,
    clean_pathway     = clean_label(feature_label),
    biological_program = assign_category(feature_label),
    treatment_direction = case_when(
      treatment_NES > 0 ~ "Up after treatment",
      treatment_NES < 0 ~ "Down after treatment",
      TRUE ~ "Zero / missing"
    ),
    opposite_disease_treatment_direction =
      !is.na(disease_NES) &
      !is.na(treatment_NES) &
      disease_NES * treatment_NES < 0
  )



all_celltypes <- sort(unique(joined$celltype))
all_programs <- sort(unique(joined$biological_program))

# =============================================================================
# UI
# =============================================================================

ui <- fluidPage(
  titlePanel("PBMC Treatment-Associated Pathway Explorer"),

  sidebarLayout(
    sidebarPanel(
      width = 3,

      h4("Recommended analysis"),

      helpText(
        "Use the three steps in order. Step 1 removes pathways consistent with recovery. ",
        "Step 2 defines the statistical treatment-not-disease candidate set. ",
        "Step 3 uses NES to prioritize the strongest remaining candidates."
      ),

      fluidRow(
        column(4, actionButton("preset_step1", "Step 1 — Recovery pre-filter", width = "100%")),
        column(4, actionButton("preset_step2", "Step 2 — FDR candidate set", width = "100%")),
        column(4, actionButton("preset_step3", "Step 3 — NES prioritization", width = "100%"))
      ),

      br(),

      uiOutput("step_guidance"),

      hr(),

      h4("Primary filters"),

      sliderInput(
        "treat_fdr",
        "Maximum treatment FDR",
        min = 0.001, max = 0.10, value = 0.05, step = 0.001
      ),

      sliderInput(
        "disease_fdr_min",
        "Minimum disease FDR",
        min = 0, max = 1, value = 0.05, step = 0.01
      ),

      checkboxInput(
        "use_disease_nes",
        "Step 3: also limit |Disease NES|",
        value = FALSE
      ),

      conditionalPanel(
        condition = "input.use_disease_nes == true",
        sliderInput(
          "disease_nes_max",
          "Maximum |Disease NES|",
          min = 0, max = 3, value = 0.75, step = 0.05
        )
      ),

      sliderInput(
        "treatment_nes_min",
        "Step 3: minimum |Treatment NES|",
        min = 0, max = 4, value = 0, step = 0.05
      ),

      hr(),

      checkboxGroupInput(
        "db",
        "Databases",
        choices = databases,
        selected = databases
      ),

      selectizeInput(
        "celltypes",
        "Cell types",
        choices = all_celltypes,
        selected = all_celltypes,
        multiple = TRUE
      ),

      selectizeInput(
        "programs",
        "Biological programs",
        choices = all_programs,
        selected = all_programs,
        multiple = TRUE
      ),

      radioButtons(
        "direction",
        "Treatment direction",
        choices = c("Both", "Up after treatment", "Down after treatment"),
        selected = "Both"
      ),

      hr(),

      checkboxInput(
        "exclude_recovery",
        "Remove pathways consistent with disease recovery",
        value = TRUE
      ),

      conditionalPanel(
        condition = "input.exclude_recovery == true",
        helpText(
          "Recovery rule: Treatment FDR is significant AND Disease NES and Treatment NES have opposite signs. ",
          "Disease FDR is intentionally not required. Thus disease-high IFN-gamma with negative treatment NES ",
          "is handled as recovery even when the PBMC untreated-vs-healthy comparison is underpowered."
        )
      ),

      hr(),

      actionButton("reset", "Reset filters", width = "100%")
    ),

    mainPanel(
      width = 9,

      fluidRow(
        column(
          3,
          wellPanel(
            h4(textOutput("n_candidates")),
            div("pathway rows in current step")
          )
        ),
        column(
          3,
          wellPanel(
            h4(textOutput("n_celltypes")),
            div("PBMC cell types")
          )
        ),
        column(
          3,
          wellPanel(
            h4(textOutput("n_programs")),
            div("biological programs")
          )
        ),
        column(
          3,
          wellPanel(
            h4(textOutput("median_nes")),
            div("median |treatment NES|")
          )
        )
      ),

      tabsetPanel(
        tabPanel(
          "Counts",
          br(),
          plotOutput("count_plot", height = "750px"),
          br(),
          DTOutput("count_table")
        ),

        tabPanel(
          "Pathways",
          br(),
          DTOutput("pathway_table")
        ),

        tabPanel(
          "Analysis funnel",
          br(),
          wellPanel(
            tags$b("Filtering order"),
            tags$ol(
              tags$li("Start: Treatment FDR < cutoff."),
              tags$li("STEP 1: remove recovery-direction pathways using opposite Disease/Treatment NES signs."),
              tags$li("STEP 2: from Step 1 survivors, require Disease FDR >= cutoff."),
              tags$li("STEP 3: from Step 2 candidates, require minimum |Treatment NES| and optionally maximum |Disease NES|."),
              tags$li("Interpret the strongest Step 3 biological programs across pathway databases.")
            ),
            tags$p(
              "This ordering maps directly onto the three analysis questions: handle expected recovery first, ",
              "define the treatment-not-disease candidate set second, then prioritize the strongest effects third."
            )
          ),
          plotOutput("funnel_plot", height = "430px"),
          DTOutput("funnel_table")
        ),

        tabPanel(
          "Recovery direction",
          br(),
          wellPanel(
            tags$b("Purpose"),
            tags$p(
              "This tab shows the pathways removed by the recovery-direction pre-filter. ",
              "Recovery is defined as Treatment FDR < the selected cutoff with Disease NES and Treatment NES ",
              "in opposite directions; Disease FDR is deliberately not required."
            ),
            tags$p(
              "Every exact pathway row listed here is removed BEFORE the treatment-not-disease screen. ",
              "Thus disease-high IFN-gamma pathways with negative treatment NES are classified as recovery even if the ",
              "PBMC untreated-vs-healthy FDR is weak. The remaining pathways are then screened for Disease FDR >= cutoff."
            )
          ),
          plotOutput("recovery_program_plot", height = "650px"),
          br(),
          DTOutput("recovery_table")
        ),

        tabPanel(
          "Step contributions",
          br(),
          wellPanel(
            tags$b("How pathway rows move through Step 1, Step 2, and Step 3"),
            tags$p(
              "This view partitions the treatment-significant pathway universe into mutually exclusive blocks ",
              "so the contribution of each analysis stage can be seen directly by cell type and database."
            ),
            tags$ul(
              tags$li(tags$b("Recovery removed before Step 1"), " — treatment-significant rows with opposite Disease/Treatment NES directions."),
              tags$li(tags$b("Step 1 only"), " — survives recovery removal but does not pass the Step 2 Disease-FDR criterion."),
              tags$li(tags$b("Step 2 only"), " — passes Step 2 but does not pass the Step 3 NES prioritization."),
              tags$li(tags$b("Step 3 retained"), " — passes recovery removal, the Step 2 FDR screen, and the Step 3 NES prioritization.")
            ),
            tags$p(
              "Because the steps are nested, these blocks are mutually exclusive. Their sum equals the treatment-significant starting universe under the current selections."
            )
          ),
          plotOutput("step_contribution_plot", height = "900px"),
          br(),
          h4("Step contributions by cell type and database"),
          DTOutput("step_contribution_table"),
          br(),
          h4("Overall contribution of each step block"),
          DTOutput("step_contribution_summary")
        ),

        tabPanel(
          "Disease vs treatment",
          br(),
          fluidRow(
            column(
              4,
              sliderInput(
                "n_scatter_labels",
                "Labels per quadrant per cell type",
                min = 0,
                max = 10,
                value = 3,
                step = 1
              )
            )
          ),
          plotOutput("scatter_plot", height = "850px")
        ),

        tabPanel(
          "Recovery audit",
          br(),
          wellPanel(
            tags$b("Sanity check"),
            tags$p(
              "This tab verifies that the exact pathway rows displayed in the Recovery Direction tab ",
              "have been removed from the downstream recovery-cleaned dataset."
            ),
            tags$p(
              "Overlap should be ZERO when recovery removal is enabled."
            )
          ),
          verbatimTextOutput("recovery_overlap_text"),
          plotOutput("recovery_program_audit_plot", height = "520px"),
          br(),
          DTOutput("recovery_program_audit_table")
        ),

        tabPanel(
          "Biological programs",
          br(),
          wellPanel(
            tags$b("Important: interpret the retained subset in the context of all strong treatment-associated pathways."),
            tags$p(
              "Plot A shows only pathways that survive the CURRENT STEP. ",
              "Step 1 = recovery-cleaned treatment pathways; Step 2 = FDR-defined candidates; ",
              "Step 3 = NES-prioritized candidates."
            ),
            tags$p(
              "Plot B starts from ALL treatment-significant pathways (Treatment FDR < 0.05). ",
              "It then shows the percentage that are also significant in untreated APECED vs healthy ",
              "versus not significant in the disease comparison."
            ),
            tags$p(
              "NES is deliberately reserved for Step 3 and does not define the Step 2 statistical candidate set."
            )
          ),
          h4("A. Programs retained by the current treatment-not-disease filter"),
          helpText(
            "IMPORTANT: Plot A is split by CELL TYPE. Point size is the number of retained pathway rows ",
            "within one cell type and biological program. Fill shows treatment direction only. ",
            "It does NOT use the blue/red disease-significance colors from Plot B."
          ),
          plotOutput("program_plot", height = "800px"),
          br(),
          h4("B. Among treatment-significant pathways, how many are also significant in disease?"),
          helpText(
            "Plot B is aggregated ACROSS ALL SELECTED CELL TYPES. It uses disease-significance classes, ",
            "not treatment direction. Therefore a large program-level percentage in Plot B can be distributed ",
            "across many small cell-type points in Plot A."
          ),
          plotOutput("program_flow_plot", height = "520px"),
          br(),
          h4("C. Program-level QC for Plot B"),
          DTOutput("program_qc_table"),
          br(),
          h4("D. Cell-type/program counts"),
          DTOutput("program_table")
        ),

        tabPanel(
          "Analysis workflow",
          br(),
          wellPanel(
            tags$b("PBMC treatment-associated pathway filtering workflow"),
            tags$p(
              "This diagram summarizes the GSEA filtering strategy and shows where ",
              "Questions 1, 2, and 3 are answered."
            )
          ),
          DiagrammeR::grVizOutput("analysis_flow_graph", width = "100%", height = "1250px")
        ),

        tabPanel(
          "IFN classification QC",
          br(),
          wellPanel(
            tags$b("Interferon pathway classification audit"),
            tags$p(
              "This table shows every interferon-labeled pathway in the selected cell types and databases, ",
              "its manually assigned biological program, and where that exact pathway row lands in the ",
              "recovery / Step 1 / Step 2 / Step 3 workflow."
            ),
            tags$p(
              "Generic Reactome parent terms such as 'Interferon signaling' should appear as ",
              tags$b("Broad / shared interferon signaling"),
              ", while explicit alpha/beta and gamma terms should remain in the type I and type II categories."
            )
          ),
          DTOutput("ifn_classification_qc_table")
        ),

        tabPanel(
          "HTML report",
          br(),
          wellPanel(
            tags$b("Download a reproducible HTML report"),
            tags$p(
              "The report uses the exact current app settings and saves the Step 1 recovery set, ",
              "Step 1 survivors, Step 2 FDR-defined candidates, Step 3 NES-prioritized candidates, ",
              "and the analysis funnel. GSEA is not rerun."
            ),
            tags$p(
              "The report explicitly maps Step 2 to Question 1, Step 3 to Question 2, and the final biological interpretation to Question 3. It also includes a disease-significant recovery-block plot showing recovered, same-direction, non-significant-treatment, and missing/not-tested contributions by cell type and database."
            ),
            downloadButton("download_html_report", "Download HTML report")
          )
        ),

        tabPanel(
          "Filter summary",
          br(),
          verbatimTextOutput("filter_text")
        )
      )
    )
  )
)

# =============================================================================
# SERVER
# =============================================================================


# Stable key used to guarantee that rows identified as recovery cannot re-enter
# downstream analyses.
joined <- joined %>%
  mutate(
    pathway_row_key = paste(celltype, database, feature_id, sep = "|||")
  )

server <- function(input, output, session) {

  analysis_stage <- reactiveVal("step2")

  output$analysis_flow_graph <- DiagrammeR::renderGrViz({
    DiagrammeR::grViz(analysis_flow_gv)
  })

  # Standard DataTable formatting used across every tab.
  standard_dt <- function(x, page_length = 25, filter_top = TRUE) {
    datatable(
      x,
      rownames = FALSE,
      filter = if (filter_top) "top" else "none",
      extensions = "Buttons",
      options = list(
        pageLength = page_length,
        scrollX = TRUE,
        autoWidth = TRUE,
        dom = "Bfrtip",
        buttons = c("copy", "csv", "excel"),
        lengthMenu = c(10, 25, 50, 100)
      )
    )
  }


  observeEvent(input$preset_step1, {
    analysis_stage("step1")
    updateSliderInput(session, "treat_fdr", value = 0.05)
    updateSliderInput(session, "disease_fdr_min", value = 0.05)
    updateCheckboxInput(session, "use_disease_nes", value = FALSE)
    updateSliderInput(session, "disease_nes_max", value = 1.00)
    updateSliderInput(session, "treatment_nes_min", value = 0)
    updateCheckboxInput(session, "exclude_recovery", value = TRUE)
  })

  observeEvent(input$preset_step2, {
    analysis_stage("step2")
    updateSliderInput(session, "treat_fdr", value = 0.05)
    updateSliderInput(session, "disease_fdr_min", value = 0.05)
    updateCheckboxInput(session, "use_disease_nes", value = FALSE)
    updateSliderInput(session, "disease_nes_max", value = 1.00)
    updateSliderInput(session, "treatment_nes_min", value = 0)
    updateCheckboxInput(session, "exclude_recovery", value = TRUE)
  })

  observeEvent(input$preset_step3, {
    analysis_stage("step3")
    updateSliderInput(session, "treat_fdr", value = 0.05)
    updateSliderInput(session, "disease_fdr_min", value = 0.05)
    updateCheckboxInput(session, "exclude_recovery", value = TRUE)
    updateSliderInput(session, "treatment_nes_min", value = 1.50)
    updateCheckboxInput(session, "use_disease_nes", value = TRUE)
    updateSliderInput(session, "disease_nes_max", value = 1.00)
  })


  observeEvent(input$reset, {
    analysis_stage("step2")
    updateSliderInput(session, "treat_fdr", value = 0.05)
    updateSliderInput(session, "disease_fdr_min", value = 0.05)
    updateCheckboxInput(session, "use_disease_nes", value = FALSE)
    updateSliderInput(session, "disease_nes_max", value = 0.75)
    updateSliderInput(session, "treatment_nes_min", value = 0)
    updateCheckboxGroupInput(session, "db", selected = databases)
    updateSelectizeInput(session, "celltypes", selected = all_celltypes)
    updateSelectizeInput(session, "programs", selected = all_programs)
    updateRadioButtons(session, "direction", selected = "Both")
    updateCheckboxInput(session, "exclude_recovery", value = TRUE)
    updateSliderInput(session, "n_scatter_labels", value = 3)
  })

  output$step_guidance <- renderUI({
    stage <- analysis_stage()

    if (stage == "step1") {
      wellPanel(
        tags$b("STEP 1 — Recovery pre-filter"),
        tags$p(
          tags$b("Purpose: "),
          "remove treatment-responsive pathways whose direction is consistent with correcting disease biology."
        ),
        tags$p(
          "Start with Treatment FDR < 0.05. Remove pathways when Disease NES and Treatment NES have opposite signs. ",
          "Disease FDR is NOT required for this recovery classification."
        ),
        tags$p(
          tags$b("Biological calibration: "),
          "disease-high IFN-gamma pathways that decrease after treatment should be captured here even when the PBMC disease comparison is weak."
        ),
        tags$p(tags$b("Step 1 output: "), "recovery-cleaned treatment-responsive pathways.")
      )
    } else if (stage == "step2") {
      wellPanel(
        tags$b("STEP 2 — FDR-defined candidate set"),
        tags$p(
          tags$b("Question 1: "),
          "from Step 1 survivors, retain pathways with Treatment FDR < 0.05 and Disease FDR >= 0.05."
        ),
        tags$p(
          "NES magnitude does NOT define membership in Step 2. This is the transparent statistical treatment-not-disease set."
        ),
        tags$p(
          tags$b("Step 2 output: "),
          "non-recovery treatment-associated pathways that are not statistically significant in the disease comparison."
        )
      )
    } else {
      wellPanel(
        tags$b("STEP 3 — NES prioritization"),
        tags$p(
          tags$b("Question 2: "),
          "from the Step 2 candidate set, prioritize larger |Treatment NES| and smaller |Disease NES|."
        ),
        tags$p(
          "Default prioritization: |Treatment NES| >= 1.5 and |Disease NES| <= 1.0. Both remain interactive."
        ),
        tags$p(
          tags$b("Question 3: "),
          "interpret the strongest Step 3 pathways/programs as candidate broader pharmacologic effects of JAK/STAT inhibition, not established side effects."
        )
      )
    }
  })

  # Main analysis funnel:
  #   Stage 1: identify treatment-significant pathways.
  #   Stage 2: identify pathways consistent with biological recovery using DIRECTION,
  #            regardless of Disease FDR:
  #              Treatment FDR < cutoff AND Disease NES * Treatment NES < 0.
  #   Stage 3: remove those recovery-direction pathways (when enabled).
  #   Stage 4: among the remaining treatment-significant pathways, retain pathways
  #            not significant in untreated APECED vs healthy.
  #
  # This order is intentional. A known disease-recovery pathway such as disease-high
  # IFN-gamma with negative Treatment NES is handled as recovery even if the small PBMC
  # untreated-vs-healthy comparison does not reach FDR significance.
  filtered <- reactive({
    x <- recovery_cleaned_treatment_universe() %>%
      mutate(
        treatment_significant = TRUE,
        disease_nonsignificant =
          !is.na(disease_FDR) &
          disease_FDR >= input$disease_fdr_min
      )

    if (analysis_stage() %in% c("step2", "step3")) {
      x <- x %>% filter(disease_nonsignificant)
    }

    if (analysis_stage() == "step3") {
      x <- x %>%
        filter(
          !is.na(abs_treatment_NES),
          abs_treatment_NES >= input$treatment_nes_min
        )

      if (isTRUE(input$use_disease_nes)) {
        x <- x %>%
          filter(
            !is.na(abs_disease_NES),
            abs_disease_NES <= input$disease_nes_max
          )
      }
    }

    x %>%
      mutate(
        analysis_stage = analysis_stage(),
        interpretation = case_when(
          analysis_stage() == "step1" ~
            "Step 1 survivor: treatment-responsive after recovery removal",
          analysis_stage() == "step2" ~
            "Step 2 candidate: recovery-cleaned + Treatment FDR significant + Disease FDR non-significant",
          analysis_stage() == "step3" ~
            "Step 3 prioritized candidate: Step 2 candidate passing NES filters",
          TRUE ~ "Filtered pathway"
        )
      ) %>%
      arrange(treatment_FDR, desc(abs_treatment_NES), abs_disease_NES)
  })

  # -----------------------------------------------------------------------
  # STEP 1 / STEP 2 / STEP 3 CONTRIBUTION BLOCKS
  # -----------------------------------------------------------------------
  # These blocks are mutually exclusive and are built from the same
  # treatment-significant starting universe used by the main workflow.
  step_contribution_blocks <- reactive({
    req(length(input$db) > 0)
    req(length(input$celltypes) > 0)
    req(length(input$programs) > 0)

    # Starting universe: treatment-significant rows under current selections.
    start <- joined %>%
      filter(
        database %in% input$db,
        celltype %in% input$celltypes,
        biological_program %in% input$programs,
        !is.na(treatment_FDR),
        treatment_FDR < input$treat_fdr
      )

    if (input$direction == "Up after treatment") {
      start <- start %>% filter(treatment_NES > 0)
    } else if (input$direction == "Down after treatment") {
      start <- start %>% filter(treatment_NES < 0)
    }

    # Recovery rows are identified exactly as in Step 1.
    recovery_keys <- recovery_direction_set() %>%
      distinct(pathway_row_key) %>%
      mutate(recovery_removed = TRUE)

    # Step 2 membership: survives recovery removal + Disease FDR non-significant.
    step2_keys <- report_step2() %>%
      distinct(pathway_row_key) %>%
      mutate(in_step2 = TRUE)

    # Step 3 membership: Step 2 + NES prioritization.
    step3_keys <- report_step3() %>%
      distinct(pathway_row_key) %>%
      mutate(in_step3 = TRUE)

    start %>%
      left_join(recovery_keys, by = "pathway_row_key") %>%
      left_join(step2_keys, by = "pathway_row_key") %>%
      left_join(step3_keys, by = "pathway_row_key") %>%
      mutate(
        recovery_removed = coalesce(recovery_removed, FALSE),
        in_step2 = coalesce(in_step2, FALSE),
        in_step3 = coalesce(in_step3, FALSE),
        step_block = case_when(
          isTRUE(input$exclude_recovery) & recovery_removed ~
            "Recovery removed before Step 1",
          in_step3 ~
            "Step 3 retained",
          in_step2 ~
            "Step 2 only",
          TRUE ~
            "Step 1 only"
        ),
        step_block = factor(
          step_block,
          levels = c(
            "Recovery removed before Step 1",
            "Step 1 only",
            "Step 2 only",
            "Step 3 retained"
          )
        )
      )
  })

  output$step_contribution_plot <- renderPlot({
    x <- step_contribution_blocks()

    validate(
      need(
        nrow(x) > 0,
        "No treatment-significant pathway rows are available under the current selections."
      )
    )

    plot_df <- x %>%
      count(database, celltype, step_block, name = "n_pathways") %>%
      group_by(database, celltype) %>%
      mutate(total = sum(n_pathways)) %>%
      ungroup()

    cell_order <- plot_df %>%
      group_by(celltype) %>%
      summarise(total_all = sum(n_pathways), .groups = "drop") %>%
      arrange(total_all) %>%
      pull(celltype)

    plot_df <- plot_df %>%
      mutate(celltype = factor(celltype, levels = cell_order))

    ggplot(
      plot_df,
      aes(
        x = n_pathways,
        y = celltype,
        fill = step_block
      )
    ) +
      geom_col(width = 0.78) +
      facet_wrap(~ database, scales = "free_x", nrow = 1) +
      labs(
        title = "Contribution of Step 1, Step 2, and Step 3 to the treatment-significant pathway universe",
        subtitle = paste0(
          "Treatment FDR < ", input$treat_fdr,
          "; Step 2 requires Disease FDR >= ", input$disease_fdr_min,
          "; Step 3 requires |Treatment NES| >= ", input$treatment_nes_min,
          if (isTRUE(input$use_disease_nes)) {
            paste0(" and |Disease NES| <= ", input$disease_nes_max)
          } else {
            ""
          },
          "."
        ),
        x = "Number of treatment-significant pathway rows",
        y = NULL,
        fill = "Analysis block"
      ) +
      theme_bw(base_size = 12) +
      theme(
        legend.position = "bottom",
        strip.text = element_text(face = "bold"),
        axis.text.y = element_text(size = 9)
      )
  })

  output$step_contribution_table <- renderDT({
    x <- step_contribution_blocks() %>%
      count(database, celltype, step_block, name = "n_pathways") %>%
      group_by(database, celltype) %>%
      mutate(
        total_treatment_significant = sum(n_pathways),
        percent_of_celltype_database =
          100 * n_pathways / total_treatment_significant
      ) %>%
      ungroup() %>%
      arrange(database, desc(total_treatment_significant), celltype, step_block)

    standard_dt(x, page_length = 25)
  })

  output$step_contribution_summary <- renderDT({
    x <- step_contribution_blocks() %>%
      count(step_block, name = "n_pathways") %>%
      mutate(
        percent = 100 * n_pathways / sum(n_pathways)
      ) %>%
      arrange(step_block)

    standard_dt(x, page_length = 10, filter_top = FALSE)
  })


  # Canonical set of rows identified as recovery.
  # These exact pathway_row_key values are removed from every downstream dataset.
  recovery_direction_set <- reactive({
    req(length(input$db) > 0)
    req(length(input$celltypes) > 0)
    req(length(input$programs) > 0)

    x <- joined %>%
      filter(
        !is.na(treatment_FDR),
        treatment_FDR < input$treat_fdr,
        !is.na(disease_NES),
        !is.na(treatment_NES),
        disease_NES * treatment_NES < 0,
        database %in% input$db,
        celltype %in% input$celltypes,
        biological_program %in% input$programs
      )

    if (input$direction != "Both") {
      x <- x %>% filter(treatment_direction == input$direction)
    }

    x %>%
      mutate(
        recovery_direction = case_when(
          disease_NES > 0 & treatment_NES < 0 ~
            "Disease-high -> decreased after treatment",
          disease_NES < 0 & treatment_NES > 0 ~
            "Disease-low -> increased after treatment",
          TRUE ~ "Opposite direction"
        ),
        recovery_disease_significance = case_when(
          is.na(disease_FDR) ~ "Disease FDR missing",
          disease_FDR < input$disease_fdr_min ~ "Disease FDR significant",
          TRUE ~ "Disease FDR not significant"
        )
      ) %>%
      distinct(pathway_row_key, .keep_all = TRUE) %>%
      arrange(treatment_FDR, desc(abs_treatment_NES), desc(abs_disease_NES))
  })


  output$n_candidates <- renderText({
    format(nrow(filtered()), big.mark = ",")
  })

  output$n_celltypes <- renderText({
    n_distinct(filtered()$celltype)
  })

  output$n_programs <- renderText({
    n_distinct(filtered()$biological_program)
  })

  output$median_nes <- renderText({
    x <- filtered()
    if (nrow(x) == 0) return("NA")
    sprintf("%.2f", median(x$abs_treatment_NES, na.rm = TRUE))
  })

  output$count_plot <- renderPlot({
    x <- filtered()
    validate(need(nrow(x) > 0, "No pathways pass the current filters."))

    plot_df <- x %>%
      count(celltype, treatment_direction, name = "n_pathways") %>%
      group_by(celltype) %>%
      mutate(total = sum(n_pathways)) %>%
      ungroup() %>%
      mutate(celltype = fct_reorder(celltype, total))

    ggplot(
      plot_df,
      aes(x = n_pathways, y = celltype, fill = treatment_direction)
    ) +
      geom_col() +
      geom_text(
        aes(label = n_pathways),
        position = position_stack(vjust = 0.5),
        size = 3
      ) +
      labs(
        title = "Pathways remaining after current filters",
        subtitle = paste0(
          if (analysis_stage() == "step1") {
            paste0(
              "STEP 1: Treatment FDR < ", input$treat_fdr,
              "; recovery-direction pathways removed"
            )
          } else if (analysis_stage() == "step2") {
            paste0(
              "STEP 2: Step 1 survivors + Disease FDR >= ", input$disease_fdr_min
            )
          } else {
            paste0(
              "STEP 3: Step 2 candidates + |Treatment NES| >= ", input$treatment_nes_min,
              if (isTRUE(input$use_disease_nes))
                paste0("; |Disease NES| <= ", input$disease_nes_max)
              else ""
            )
          }
        ),
        x = "Number of pathway rows",
        y = NULL,
        fill = NULL
      ) +
      theme_bw(base_size = 12) +
      theme(legend.position = "top")
  })

  output$count_table <- renderDT({
    x <- filtered() %>%
      count(celltype, database, treatment_direction, name = "n_pathways") %>%
      arrange(desc(n_pathways))

    standard_dt(x, page_length = 25)
  })

  output$pathway_table <- renderDT({
    x <- filtered() %>%
      select(
        celltype, database, feature_id, clean_pathway, biological_program,
        treatment_direction,
        treatment_NES, treatment_FDR,
        disease_NES, abs_disease_NES, disease_FDR,
        analysis_stage,
        treatment_significant,
        opposite_disease_treatment_direction,
        disease_nonsignificant,
        interpretation
      )

    x <- x %>%
      rename(
        recovery_direction_flag = opposite_disease_treatment_direction
      )

    standard_dt(x, page_length = 25)
  })

  output$scatter_plot <- renderPlot({
    x <- filtered()
    validate(need(nrow(x) > 0, "No pathways pass the current filters."))

    # Define quadrant from the signs of disease NES (x) and treatment NES (y).
    # Labels are selected independently WITHIN EACH CELL TYPE AND QUADRANT.
    x <- x %>%
      mutate(
        quadrant = case_when(
          disease_NES >= 0 & treatment_NES >= 0 ~ "Q1: disease + / treatment +",
          disease_NES <  0 & treatment_NES >= 0 ~ "Q2: disease - / treatment +",
          disease_NES <  0 & treatment_NES <  0 ~ "Q3: disease - / treatment -",
          disease_NES >= 0 & treatment_NES <  0 ~ "Q4: disease + / treatment -",
          TRUE ~ NA_character_
        )
      )

    n_per_quad <- input$n_scatter_labels

    label_df <- x %>%
      filter(!is.na(quadrant), !is.na(clean_pathway)) %>%
      group_by(celltype, quadrant) %>%
      arrange(
        treatment_FDR,
        desc(abs_treatment_NES),
        abs_disease_NES,
        .by_group = TRUE
      ) %>%
      slice_head(n = n_per_quad) %>%
      ungroup()

    p <- ggplot(x, aes(x = disease_NES, y = treatment_NES)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
      geom_point(aes(shape = database), alpha = 0.70, size = 2) +
      facet_wrap(~ celltype, ncol = 4, scales = "free") +
      labs(
        title = "Disease NES vs paired treatment NES",
        subtitle = paste0(
          n_per_quad,
          " label", ifelse(n_per_quad == 1, "", "s"),
          " per quadrant per cell type; ranked by treatment FDR then |treatment NES|"
        ),
        x = "Mild untreated vs healthy — NES",
        y = "Mild treated vs untreated — NES",
        shape = "Database"
      ) +
      theme_bw(base_size = 11) +
      theme(
        legend.position = "bottom",
        strip.text = element_text(face = "bold", size = 8)
      )

    if (n_per_quad > 0 && nrow(label_df) > 0) {
      p <- p +
        ggrepel::geom_text_repel(
          data = label_df,
          aes(label = clean_pathway),
          size = 2.4,
          max.overlaps = Inf,
          box.padding = 0.25,
          point.padding = 0.15,
          min.segment.length = 0,
          seed = 42,
          show.legend = FALSE
        )
    }

    p
  })

  # -----------------------------------------------------------------------
  # CANONICAL RECOVERY-CLEANED UNIVERSE
  # -----------------------------------------------------------------------
  # All downstream analysis tabs must start here.
  #
  # 1. Treatment FDR < selected cutoff
  # 2. Identify recovery by opposite Disease/Treatment NES signs
  # 3. Remove those recovery-direction pathways when recovery removal is ON
  #
  # Disease FDR is NOT used to construct this universe.
  recovery_cleaned_treatment_universe <- reactive({
    req(length(input$db) > 0)
    req(length(input$celltypes) > 0)
    req(length(input$programs) > 0)

    x <- joined %>%
      filter(
        !is.na(treatment_FDR),
        treatment_FDR < input$treat_fdr,
        database %in% input$db,
        celltype %in% input$celltypes,
        biological_program %in% input$programs
      )

    if (input$direction != "Both") {
      x <- x %>% filter(treatment_direction == input$direction)
    }

    # IMPORTANT: remove the EXACT rows shown in the Recovery Direction tab.
    # This is an anti-join on a stable key, not a second re-evaluation of the rule.
    if (isTRUE(input$exclude_recovery)) {
      recovery_keys <- recovery_direction_set() %>%
        distinct(pathway_row_key)

      x <- x %>%
        anti_join(recovery_keys, by = "pathway_row_key")
    }

    x
  })



  # Program universe used as the denominator for biological-program summaries.
  # IMPORTANT: this is AFTER recovery removal.
  # Treatment NES threshold is retained only where the user explicitly requests it.
  treatment_program_universe <- reactive({
    x <- recovery_cleaned_treatment_universe()

    if (analysis_stage() == "step3") {
      x <- x %>%
        filter(
          !is.na(abs_treatment_NES),
          abs_treatment_NES >= input$treatment_nes_min
        )

      if (isTRUE(input$use_disease_nes)) {
        x <- x %>%
          filter(
            !is.na(abs_disease_NES),
            abs_disease_NES <= input$disease_nes_max
          )
      }
    }

    x
  })


  # Recovery-cleaned treatment-significant universe.
  # No Treatment NES magnitude threshold.
  treatment_significant_program_universe <- reactive({
    recovery_cleaned_treatment_universe()
  })


  # Disease-FDR context among pathways that SURVIVED recovery removal.
  treatment_significant_flow_rows <- reactive({
    x <- recovery_cleaned_treatment_universe()

    if (nrow(x) == 0) return(tibble())

    x %>%
      mutate(
        disease_significance_class = case_when(
          is.na(disease_FDR) ~ "Missing disease comparison",
          disease_FDR < input$disease_fdr_min ~ "Also significant in disease",
          disease_FDR >= input$disease_fdr_min ~ "Not significant in disease",
          TRUE ~ "Missing disease comparison"
        )
      )
  })


  recovery_overlap_audit <- reactive({
    removed <- recovery_direction_set() %>%
      distinct(pathway_row_key, biological_program, celltype, database, feature_id)

    surviving <- recovery_cleaned_treatment_universe() %>%
      distinct(pathway_row_key, biological_program, celltype, database, feature_id)

    overlap <- inner_join(
      removed,
      surviving,
      by = "pathway_row_key",
      suffix = c("_removed", "_surviving")
    )

    list(
      removed = removed,
      surviving = surviving,
      overlap = overlap
    )
  })

  output$recovery_overlap_text <- renderText({
    a <- recovery_overlap_audit()

    paste0(
      "Recovery rows identified: ", nrow(a$removed), "\n",
      "Treatment-significant rows after recovery removal: ", nrow(a$surviving), "\n",
      "Exact recovery rows still present downstream: ", nrow(a$overlap), "\n\n",
      if (isTRUE(input$exclude_recovery) && nrow(a$overlap) == 0) {
        "PASS: every row shown in the Recovery Direction tab has been removed."
      } else if (!isTRUE(input$exclude_recovery)) {
        "Recovery removal is currently disabled."
      } else {
        "FAIL: some recovery rows are still present downstream."
      }
    )
  })

  recovery_program_audit <- reactive({
    treatment_total <- joined %>%
      filter(
        !is.na(treatment_FDR),
        treatment_FDR < input$treat_fdr,
        database %in% input$db,
        celltype %in% input$celltypes,
        biological_program %in% input$programs
      )

    if (input$direction != "Both") {
      treatment_total <- treatment_total %>%
        filter(treatment_direction == input$direction)
    }

    total_by_program <- treatment_total %>%
      filter(biological_program != "Other") %>%
      count(biological_program, name = "treatment_significant_total")

    removed_by_program <- recovery_direction_set() %>%
      filter(biological_program != "Other") %>%
      count(biological_program, name = "recovery_removed")

    surviving_by_program <- recovery_cleaned_treatment_universe() %>%
      filter(biological_program != "Other") %>%
      count(biological_program, name = "surviving_after_recovery")

    total_by_program %>%
      full_join(removed_by_program, by = "biological_program") %>%
      full_join(surviving_by_program, by = "biological_program") %>%
      mutate(
        across(
          c(treatment_significant_total, recovery_removed, surviving_after_recovery),
          ~coalesce(.x, 0L)
        ),
        pct_removed_as_recovery =
          100 * recovery_removed / pmax(treatment_significant_total, 1L)
      ) %>%
      arrange(desc(recovery_removed))
  })

  output$recovery_program_audit_plot <- renderPlot({
    x <- recovery_program_audit() %>%
      select(
        biological_program,
        recovery_removed,
        surviving_after_recovery
      ) %>%
      pivot_longer(
        cols = c(recovery_removed, surviving_after_recovery),
        names_to = "status",
        values_to = "n_pathways"
      ) %>%
      mutate(
        status = recode(
          status,
          recovery_removed = "Removed as recovery",
          surviving_after_recovery = "Surviving after recovery removal"
        )
      )

    validate(need(nrow(x) > 0, "No pathways available for audit."))

    ggplot(
      x,
      aes(
        x = n_pathways,
        y = fct_reorder(biological_program, n_pathways, .fun = sum),
        fill = status
      )
    ) +
      geom_col() +
      labs(
        title = "What was removed vs what actually survives recovery filtering",
        subtitle = "A biological program can still appear downstream if it contains other non-recovery pathway rows.",
        x = "Number of pathway rows",
        y = NULL,
        fill = NULL
      ) +
      theme_bw(base_size = 11) +
      theme(legend.position = "bottom")
  })

  output$recovery_program_audit_table <- renderDT({
    standard_dt(recovery_program_audit(), page_length = 25)
  })

  analysis_funnel <- reactive({
    base <- joined %>%
      filter(
        database %in% input$db,
        celltype %in% input$celltypes,
        biological_program %in% input$programs
      )

    if (input$direction != "Both") {
      base <- base %>% filter(treatment_direction == input$direction)
    }

    treatment_sig <- base %>%
      filter(!is.na(treatment_FDR), treatment_FDR < input$treat_fdr)

    recovery <- recovery_direction_set()
    step1_survivors <- recovery_cleaned_treatment_universe()

    step2_candidates <- step1_survivors %>%
      filter(!is.na(disease_FDR), disease_FDR >= input$disease_fdr_min)

    step3_candidates <- step2_candidates %>%
      filter(!is.na(abs_treatment_NES), abs_treatment_NES >= input$treatment_nes_min)

    if (isTRUE(input$use_disease_nes)) {
      step3_candidates <- step3_candidates %>%
        filter(!is.na(abs_disease_NES), abs_disease_NES <= input$disease_nes_max)
    }

    current <- switch(
      analysis_stage(),
      step1 = step1_survivors,
      step2 = step2_candidates,
      step3 = step3_candidates,
      step2_candidates
    )

    tibble(
      stage = c(
        "Treatment FDR significant",
        "Recovery rows removed",
        "STEP 1 survivors",
        "STEP 2: + Disease FDR non-significant",
        "STEP 3: + NES prioritization",
        "Current displayed result"
      ),
      n_pathway_rows = c(
        nrow(treatment_sig),
        nrow(recovery),
        nrow(step1_survivors),
        nrow(step2_candidates),
        nrow(step3_candidates),
        nrow(current)
      )
    )
  })

  output$funnel_plot <- renderPlot({
    x <- analysis_funnel()

    ggplot(
      x,
      aes(
        x = n_pathway_rows,
        y = forcats::fct_rev(factor(stage, levels = stage))
      )
    ) +
      geom_col() +
      geom_text(
        aes(label = n_pathway_rows),
        hjust = -0.15,
        size = 4
      ) +
      scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
      labs(
        title = "Analysis funnel",
        subtitle = "Step 1 removes recovery; Step 2 applies the FDR candidate rule; Step 3 applies NES prioritization.",
        x = "Pathway rows",
        y = NULL
      ) +
      theme_bw(base_size = 11)
  })

  output$funnel_table <- renderDT({
    standard_dt(analysis_funnel(), page_length = 25)
  })

  output$recovery_program_plot <- renderPlot({
    x <- recovery_direction_set() %>%
      filter(biological_program != "Other") %>%
      count(biological_program, recovery_direction, name = "n_pathways")

    validate(need(nrow(x) > 0, "No opposite-direction treatment-significant pathways for current selections."))

    ggplot(
      x,
      aes(
        x = n_pathways,
        y = fct_reorder(biological_program, n_pathways, .fun = sum),
        fill = recovery_direction
      )
    ) +
      geom_col() +
      labs(
        title = "Recovery-direction pathways identified before the treatment-not-disease screen",
        subtitle = paste0(
          "Treatment FDR < ", input$treat_fdr,
          "; Disease FDR is not used in this recovery-direction view."
        ),
        x = "Number of pathway rows",
        y = NULL,
        fill = "Recovery direction"
      ) +
      theme_bw(base_size = 11) +
      theme(legend.position = "bottom")
  })

  output$recovery_table <- renderDT({
    x <- recovery_direction_set() %>%
      select(
        celltype,
        database,
        feature_id,
        clean_pathway,
        biological_program,
        recovery_direction,
        recovery_disease_significance,
        treatment_NES,
        treatment_FDR,
        disease_NES,
        disease_FDR
      )

    standard_dt(x, page_length = 25)
  })

  program_summary <- reactive({
    retained <- filtered()
    universe <- treatment_program_universe()

    if (nrow(universe) == 0) return(tibble())

    universe_sum <- universe %>%
      filter(biological_program != "Other") %>%
      group_by(celltype, biological_program) %>%
      summarise(
        n_recovery_cleaned_treatment_associated = n(),
        .groups = "drop"
      )

    retained_sum <- retained %>%
      filter(biological_program != "Other") %>%
      group_by(celltype, biological_program) %>%
      summarise(
        n_retained = n(),
        n_up = sum(treatment_NES > 0, na.rm = TRUE),
        n_down = sum(treatment_NES < 0, na.rm = TRUE),
        n_databases = n_distinct(database),
        n_GO_BP = sum(database == "GO_BP"),
        n_Reactome = sum(database == "Reactome"),
        n_Hallmark = sum(database == "Hallmark"),
        median_treatment_NES = median(treatment_NES, na.rm = TRUE),
        median_abs_treatment_NES = median(abs_treatment_NES, na.rm = TRUE),
        median_disease_NES = median(disease_NES, na.rm = TRUE),
        median_abs_disease_NES = median(abs_disease_NES, na.rm = TRUE),
        best_treatment_FDR = min(treatment_FDR, na.rm = TRUE),
        .groups = "drop"
      )

    universe_sum %>%
      left_join(
        retained_sum,
        by = c("celltype", "biological_program")
      ) %>%
      mutate(
        n_retained = coalesce(n_retained, 0L),
        n_up = coalesce(n_up, 0L),
        n_down = coalesce(n_down, 0L),
        n_databases = coalesce(n_databases, 0L),
        n_GO_BP = coalesce(n_GO_BP, 0L),
        n_Reactome = coalesce(n_Reactome, 0L),
        n_Hallmark = coalesce(n_Hallmark, 0L),
        retained_pct_of_treatment =
          100 * n_retained / pmax(n_recovery_cleaned_treatment_associated, 1L),
        dominant_treatment_direction = case_when(
          n_retained == 0 ~ "No retained pathways",
          n_up > n_down ~ "Mostly up after treatment",
          n_down > n_up ~ "Mostly down after treatment",
          n_up == n_down & n_up > 0 ~ "Mixed / tied",
          TRUE ~ "No direction"
        ),
        multi_database = n_databases >= 2
      ) %>%
      arrange(
        desc(n_retained),
        desc(n_databases),
        celltype,
        biological_program
      )
  })

  output$program_plot <- renderPlot({
    x <- program_summary() %>%
      filter(n_retained > 0)

    validate(need(nrow(x) > 0, "No biological programs pass the current filters."))

    # Order biological programs by how broadly they occur across PBMC cell types.
    breadth <- x %>%
      group_by(biological_program) %>%
      summarise(
        n_celltypes = n_distinct(celltype),
        total_retained = sum(n_retained),
        .groups = "drop"
      ) %>%
      arrange(desc(n_celltypes), desc(total_retained)) %>%
      mutate(
        program_label = paste0(
          biological_program, " [", n_celltypes, "/",
          n_distinct(x$celltype), " cell types]"
        )
      )

    cell_order <- x %>%
      group_by(celltype) %>%
      summarise(total = sum(n_retained), .groups = "drop") %>%
      arrange(total) %>%
      pull(celltype)

    x <- x %>%
      left_join(
        breadth %>% select(biological_program, program_label),
        by = "biological_program"
      ) %>%
      mutate(
        celltype = factor(celltype, levels = cell_order),
        program_label = factor(
          program_label,
          levels = rev(breadth$program_label)
        ),
        dominant_treatment_direction = factor(
          dominant_treatment_direction,
          levels = c(
            "Mostly down after treatment",
            "Mixed / tied",
            "Mostly up after treatment"
          )
        )
      )

    ggplot(
      x,
      aes(
        x = program_label,
        y = celltype,
        size = n_retained,
        fill = dominant_treatment_direction,
        shape = multi_database
      )
    ) +
      geom_point(alpha = 0.9, stroke = 0.9, colour = "grey20") +
      scale_size_area(max_size = 10) +
      scale_shape_manual(values = c(`FALSE` = 21, `TRUE` = 24)) +
      labs(
        title = "Biological programs retained after recovery removal and current-step filtering",
        subtitle = paste0(
          "Programs are ordered by PBMC cell-type breadth; labels show how many selected cell types contain the program. ",
          "Size = retained pathway rows (secondary support); fill = direction; triangle = >=2 databases."
        ),
        x = NULL,
        y = NULL,
        size = "Retained pathway rows",
        fill = "Direction among retained rows",
        shape = "Program represented in >=2 databases"
      ) +
      theme_bw(base_size = 11) +
      theme(
        axis.text.x = element_text(angle = 35, hjust = 1),
        legend.position = "bottom",
        panel.grid.minor = element_blank()
      )
  })

  output$program_flow_plot <- renderPlot({
    x <- treatment_significant_flow_rows() %>%
      filter(biological_program != "Other") %>%
      count(biological_program, disease_significance_class, name = "n_pathway_rows") %>%
      group_by(biological_program) %>%
      mutate(
        total_treatment_significant = sum(n_pathway_rows),
        pct = 100 * n_pathway_rows / total_treatment_significant
      ) %>%
      ungroup()

    validate(need(nrow(x) > 0, "No treatment-significant pathways for the current selections."))

    program_order <- x %>%
      group_by(biological_program) %>%
      summarise(
        total = sum(n_pathway_rows),
        pct_also_disease_sig = sum(
          pct[disease_significance_class == "Also significant in disease"],
          na.rm = TRUE
        ),
        .groups = "drop"
      ) %>%
      arrange(pct_also_disease_sig) %>%
      pull(biological_program)

    x <- x %>%
      mutate(
        biological_program = factor(biological_program, levels = program_order)
      )

    totals <- x %>%
      distinct(biological_program, total_treatment_significant)

    ggplot(
      x,
      aes(
        x = pct,
        y = biological_program,
        fill = disease_significance_class
      )
    ) +
      geom_col(width = 0.75) +
      scale_fill_manual(
        values = c(
          "Also significant in disease" = "grey35",
          "Not significant in disease" = "grey75",
          "Missing disease comparison" = "white"
        ),
        drop = FALSE
      ) +
      geom_text(
        data = totals,
        aes(
          x = 101,
          y = biological_program,
          label = paste0("n=", total_treatment_significant)
        ),
        inherit.aes = FALSE,
        hjust = 0,
        size = 3
      ) +
      coord_cartesian(xlim = c(0, 108), clip = "off") +
      labs(
        title = "Recovery-cleaned treatment-significant pathways: disease significance status",
        subtitle = paste0(
          "Aggregated across all selected cell types. Denominator: treatment-significant pathways AFTER recovery removal; ",
          "no Treatment NES cutoff. Disease significance is defined by Disease FDR < ",
          input$disease_fdr_min, "."
        ),
        x = "Percent of treatment-significant pathway rows",
        y = NULL,
        fill = "Disease comparison"
      ) +
      scale_x_continuous(
        breaks = seq(0, 100, 20),
        labels = function(z) paste0(z, "%")
      ) +
      theme_bw(base_size = 11) +
      theme(
        legend.position = "bottom",
        plot.margin = margin(5.5, 45, 5.5, 5.5)
      )
  })

  output$program_qc_table <- renderDT({
    x <- treatment_significant_flow_rows() %>%
      filter(biological_program != "Other") %>%
      group_by(biological_program) %>%
      summarise(
        recovery_cleaned_treatment_significant_total = n(),
        also_significant_in_disease =
          sum(disease_significance_class == "Also significant in disease"),
        not_significant_in_disease =
          sum(disease_significance_class == "Not significant in disease"),
        missing_disease =
          sum(disease_significance_class == "Missing disease comparison"),
        pct_also_significant_in_disease =
          100 * also_significant_in_disease / pmax(recovery_cleaned_treatment_significant_total, 1L),
        pct_not_significant_in_disease =
          100 * not_significant_in_disease / pmax(recovery_cleaned_treatment_significant_total, 1L),
        .groups = "drop"
      ) %>%
      arrange(desc(recovery_cleaned_treatment_significant_total))

    standard_dt(x, page_length = 25) %>%
      formatRound(
        columns = c(
          "pct_also_significant_in_disease",
          "pct_not_significant_in_disease"
        ),
        digits = 1
      )
  })

  output$program_table <- renderDT({
    context <- treatment_significant_flow_rows() %>%
      filter(biological_program != "Other") %>%
      group_by(celltype, biological_program) %>%
      summarise(
        recovery_cleaned_treatment_significant_n = n(),
        also_disease_significant_n =
          sum(disease_significance_class == "Also significant in disease"),
        disease_nonsignificant_n =
          sum(disease_significance_class == "Not significant in disease"),
        missing_disease_n =
          sum(disease_significance_class == "Missing disease comparison"),
        pct_also_disease_significant =
          100 * also_disease_significant_n / pmax(recovery_cleaned_treatment_significant_n, 1L),
        pct_disease_nonsignificant =
          100 * disease_nonsignificant_n / pmax(recovery_cleaned_treatment_significant_n, 1L),
        .groups = "drop"
      )

    retained_info <- program_summary() %>%
      select(
        celltype,
        biological_program,
        n_retained,
        retained_pct_of_treatment,
        n_up,
        n_down,
        dominant_treatment_direction,
        n_databases,
        n_GO_BP,
        n_Reactome,
        n_Hallmark,
        median_treatment_NES,
        median_abs_treatment_NES,
        median_disease_NES,
        median_abs_disease_NES,
        best_treatment_FDR
      )

    x <- context %>%
      left_join(retained_info, by = c("celltype", "biological_program"))

    standard_dt(x, page_length = 25) %>%
      formatRound(
        columns = c(
          "pct_also_disease_significant",
          "pct_disease_nonsignificant",
          "retained_pct_of_treatment",
          "median_treatment_NES",
          "median_abs_treatment_NES",
          "median_disease_NES",
          "median_abs_disease_NES"
        ),
        digits = 2
      ) %>%
      formatSignif(columns = "best_treatment_FDR", digits = 3)
  })


  # -----------------------------------------------------------------------
  # SELF-CONTAINED QUARTO TEMPLATE
  # -----------------------------------------------------------------------
  # v32:
  # - removes duplicated Question 1/2/3 block from top of Shiny app
  # - removes duplicated full question list from top of HTML report
  # - questions now appear only at the stages where they are answered
  # - retains small-cohort rationale and all prior analysis features
  generate_pbmc_qmd <- function() {
    c(
    "---",
    "title: \"PBMC Treatment-Associated Pathway Report\"",
    "subtitle: \"Recovery-aware three-step pathway analysis\"",
    "date: today",
    "format:",
    "  html:",
    "    toc: true",
    "    toc-depth: 3",
    "    number-sections: true",
    "    embed-resources: true",
    "    code-fold: true",
    "    code-tools: true",
    "    theme: cosmo",
    "    page-layout: full",
    "    fig-width: 16",
    "    fig-height: 8",
    "    fig-dpi: 120",
    "params:",
    "  report_rds: null",
    "execute:",
    "  echo: false",
    "  warning: false",
    "  message: false",
    "---",
    "",
    "<style>",
    ".quarto-container.page-columns {",
    "  max-width: 1700px !important;",
    "}",
    ".cell-output-display {",
    "  overflow-x: auto;",
    "}",
    ".dataTables_wrapper {",
    "  font-size: 0.92rem;",
    "}",
    ".dataTables_wrapper table.dataTable {",
    "  width: 100% !important;",
    "}",
    ".report-comment {",
    "  border-left: 4px solid #6c757d;",
    "  background: #f7f7f7;",
    "  padding: 0.8rem 1rem;",
    "  margin: 1rem 0 1.25rem 0;",
    "}",
    "</style>",
    "",
    "```{r}",
    "suppressPackageStartupMessages({",
    "  library(dplyr)",
    "  library(tidyr)",
    "  library(stringr)",
    "  library(ggplot2)",
    "  library(ggrepel)",
    "  library(forcats)",
    "  library(DT)",
    "  library(DiagrammeR)",
    "})",
    "",
    "if (is.null(params$report_rds) || !file.exists(params$report_rds)) {",
    "  stop(\"Report data RDS was not supplied or does not exist.\")",
    "}",
    "",
    "report <- readRDS(params$report_rds)",
    "",
    "settings <- report$settings",
    "recovery <- report$recovery",
    "step_contribution_blocks <- report$step_contribution_blocks",
    "step1 <- report$step1",
    "step2 <- report$step2",
    "step3 <- report$step3",
    "funnel <- report$funnel",
    "ifn_classification_qc <- report$interferon_classification_qc",
    "analysis_flow_gv <- report$analysis_flow_gv",
    "",
    "report_dt <- function(x, page_length = 25, filter_top = TRUE, round_digits = NULL) {",
    "  if (!is.null(round_digits)) {",
    "    num_cols <- names(x)[vapply(x, is.numeric, logical(1))]",
    "    for (nm in num_cols) {",
    "      x[[nm]] <- round(x[[nm]], round_digits)",
    "    }",
    "  }",
    "",
    "  DT::datatable(",
    "    x,",
    "    rownames = FALSE,",
    "    filter = if (filter_top) \"top\" else \"none\",",
    "    extensions = c(\"Buttons\"),",
    "    options = list(",
    "      pageLength = page_length,",
    "      scrollX = TRUE,",
    "      autoWidth = TRUE,",
    "      dom = \"Bfrtip\",",
    "      buttons = c(\"copy\", \"csv\", \"excel\"),",
    "      lengthMenu = c(10, 25, 50, 100)",
    "    )",
    "  )",
    "}",
    "",
    "comment_box <- function(text) {",
    "  cat(",
    "    paste0(",
    "      '<div class=\"report-comment\"><strong>Interpretation:</strong> ',",
    "      text,",
    "      '</div>'",
    "    )",
    "  )",
    "}",
    "",
    "fmt_n <- function(x) format(x, big.mark = \",\", scientific = FALSE)",
    "```",
    "",
    "# Analysis workflow",
    "",
    "The custom PBMC pathway-filtering app applies the filtering strategy summarized below.",
    "",
    "```{r}",
    "#| column: page",
    "#| out-width: 100%",
    "DiagrammeR::grViz(analysis_flow_gv, width = 1600, height = 1500)",
    "```",
    "",
    "# Workflow overview",
    "",
    "The report presents each research question only at the stage where it is addressed:",
    "",
    "- **Step 1 — Recovery pre-filter:** removes treatment responses consistent with disease recovery before candidate identification.",
    "- **Step 2 — FDR-defined candidate set:** addresses **Question 1**.",
    "- **Step 3 — NES prioritization:** addresses **Question 2**.",
    "- **Biological interpretation of Step 3:** addresses **Question 3**.",
    "",
    "# Statistical interpretation with a small PBMC cohort",
    "",
    "The untreated-versus-healthy PBMC comparison is based on only **4 subjects** and is therefore expected to have limited statistical power. For that reason, the workflow does not interpret a non-significant disease FDR as evidence that a pathway is unaffected by disease.",
    "",
    "The analysis uses the following logic:",
    "",
    "1. **Treatment response:** require Treatment FDR below the selected cutoff.",
    "2. **Remove recovery:** if Disease NES and Treatment NES are opposite in direction, classify the pathway as recovery regardless of Disease FDR.",
    "3. **Candidate non-recovery effects:** among the remaining pathways, use Disease FDR above the selected cutoff as one criterion, but interpret this only as **not statistically significant in this small PBMC comparison**.",
    "4. **Prioritize:** use large |Treatment NES| and, when enabled, small |Disease NES| to focus on the strongest effects least suggestive of underlying disease biology.",
    "5. **Biological confidence:** emphasize coherent biological programs and cross-database support rather than expecting perfect statistical separation from 4 subjects.",
    "",
    "> **Because the untreated-versus-healthy PBMC comparison is based on a small number of subjects and is underpowered, lack of disease significance was not interpreted as evidence of no disease effect. Recovery was therefore identified using disease-versus-treatment directionality, while residual disease NES magnitude was used to help prioritize non-recovery treatment-associated pathways.**",
    "",
    "# Interferon classification QC",
    "",
    "This audit makes the manual interferon grouping transparent. Generic parent terms such as **Interferon signaling** are assigned to **Broad / shared interferon signaling**; explicit IFN-alpha/beta terms remain in the type I category and explicit IFN-gamma terms remain in the type II category.",
    "",
    "The `classification_stage` column shows where each exact interferon-labeled pathway row lands in the current recovery / Step 1 / Step 2 / Step 3 workflow.",
    "",
    "```{r}",
    "if (is.null(ifn_classification_qc) || nrow(ifn_classification_qc) == 0) {",
    "  cat(\"No interferon-labeled pathways were found for the selected cell types/databases.\")",
    "} else {",
    "  report_dt(",
    "    ifn_classification_qc,",
    "    page_length = 50,",
    "    filter_top = TRUE,",
    "    round_digits = 3",
    "  )",
    "}",
    "```",
    "",
    "```{r, results='asis'}",
    "if (!is.null(ifn_classification_qc) && nrow(ifn_classification_qc) > 0) {",
    "  ifn_qc_summary <- ifn_classification_qc %>%",
    "    count(biological_program, classification_stage, name = \"n_rows\") %>%",
    "    arrange(biological_program, desc(n_rows))",
    "",
    "  broad_n <- ifn_classification_qc %>%",
    "    filter(biological_program == \"Broad / shared interferon signaling\") %>%",
    "    nrow()",
    "",
    "  comment_box(",
    "    paste0(",
    "      \"The classification audit contains \", nrow(ifn_classification_qc),",
    "      \" interferon-labeled pathway rows. \",",
    "      broad_n, \" are assigned to the Broad / shared interferon signaling category. \",",
    "      \"A zero count here means that no generic interferon parent pathways are present in the selected data/settings, rather than that the category is absent from the classifier.\"",
    "    )",
    "  )",
    "}",
    "```",
    "",
    "# Current analysis settings",
    "",
    "```{r}",
    "settings_tbl <- tibble(",
    "  Setting = c(",
    "    \"Treatment FDR cutoff\",",
    "    \"Disease FDR minimum\",",
    "    \"Recovery removal\",",
    "    \"Recovery rule\",",
    "    \"Step 3 minimum |Treatment NES|\",",
    "    \"Step 3 use maximum |Disease NES|\",",
    "    \"Step 3 maximum |Disease NES|\",",
    "    \"Databases\",",
    "    \"Treatment direction\",",
    "    \"Selected cell types\",",
    "    \"Selected biological programs\"",
    "  ),",
    "  Value = c(",
    "    as.character(settings$treat_fdr),",
    "    as.character(settings$disease_fdr_min),",
    "    as.character(settings$exclude_recovery),",
    "    \"Treatment FDR significant + opposite Disease/Treatment NES signs\",",
    "    as.character(settings$treatment_nes_min),",
    "    as.character(settings$use_disease_nes),",
    "    if (isTRUE(settings$use_disease_nes)) as.character(settings$disease_nes_max) else \"Not applied\",",
    "    paste(settings$db, collapse = \", \"),",
    "    settings$direction,",
    "    paste(settings$celltypes, collapse = \", \"),",
    "    paste(settings$programs, collapse = \", \")",
    "  )",
    ")",
    "",
    "report_dt(settings_tbl, page_length = 25, filter_top = FALSE)",
    "```",
    "",
    "# Analysis funnel",
    "",
    "```{r}",
    "#| fig-width: 13",
    "#| fig-height: 5",
    "ggplot(",
    "  funnel,",
    "  aes(",
    "    x = n_pathway_rows,",
    "    y = forcats::fct_rev(factor(stage, levels = stage))",
    "  )",
    ") +",
    "  geom_col() +",
    "  geom_text(",
    "    aes(label = fmt_n(n_pathway_rows)),",
    "    hjust = -0.15,",
    "    size = 4",
    "  ) +",
    "  scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +",
    "  labs(",
    "    title = \"Three-step analysis funnel\",",
    "    subtitle = \"Recovery removal precedes the FDR candidate screen; NES is reserved for prioritization.\",",
    "    x = \"Pathway rows\",",
    "    y = NULL",
    "  ) +",
    "  theme_bw(base_size = 12)",
    "```",
    "",
    "```{r}",
    "report_dt(funnel, page_length = 25, filter_top = FALSE)",
    "```",
    "",
    "```{r, results='asis'}",
    "funnel_lookup <- setNames(funnel$n_pathway_rows, funnel$stage)",
    "n_start <- funnel_lookup[[\"Treatment FDR significant\"]]",
    "n_recovery <- funnel_lookup[[\"Recovery rows removed\"]]",
    "n_step1 <- funnel_lookup[[\"STEP 1 survivors\"]]",
    "n_step2 <- funnel_lookup[[\"STEP 2: + Disease FDR non-significant\"]]",
    "n_step3 <- funnel_lookup[[\"STEP 3: + NES prioritization\"]]",
    "",
    "if (all(!is.null(c(n_start, n_recovery, n_step1, n_step2, n_step3)))) {",
    "  recovery_pct <- ifelse(n_start > 0, 100 * n_recovery / n_start, NA_real_)",
    "  step2_pct <- ifelse(n_step1 > 0, 100 * n_step2 / n_step1, NA_real_)",
    "  step3_pct <- ifelse(n_step2 > 0, 100 * n_step3 / n_step2, NA_real_)",
    "",
    "  comment_box(",
    "    paste0(",
    "      \"The treatment-significant starting set contains \", fmt_n(n_start), \" pathway rows. \",",
    "      fmt_n(n_recovery), \" (\", sprintf(\"%.1f\", recovery_pct), \"%) are classified as recovery by direction and removed in Step 1. \",",
    "      \"Step 2 retains \", fmt_n(n_step2), \" rows (\", sprintf(\"%.1f\", step2_pct),",
    "      \"% of Step 1 survivors) after requiring disease FDR to be non-significant. \",",
    "      \"Step 3 retains \", fmt_n(n_step3), \" rows (\", sprintf(\"%.1f\", step3_pct),",
    "      \"% of Step 2) after the selected NES prioritization thresholds.\"",
    "    )",
    "  )",
    "}",
    "```",
    "",
    "# Step 1 — Recovery pre-filter",
    "",
    "> **Purpose:** This step does not directly answer one of the three research questions. It is a required preprocessing step that prevents genuine disease-recovery biology from being misclassified as a broader treatment effect in Questions 1–3.",
    "",
    "Step 1 begins with pathways significantly altered after treatment and removes pathways whose treatment direction is opposite to the disease NES. Disease FDR is deliberately **not required** for this recovery classification.",
    "",
    "This prevents plausible disease-recovery biology from entering the downstream treatment-not-disease candidate set solely because the PBMC untreated-versus-healthy comparison is weak.",
    "",
    "## Recovery-direction pathways removed",
    "",
    "```{r}",
    "#| fig-width: 14",
    "#| fig-height: 7",
    "if (nrow(recovery) == 0) {",
    "  cat(\"No pathways met the recovery-direction rule under the selected settings.\")",
    "} else {",
    "  recovery_plot <- recovery %>%",
    "    filter(biological_program != \"Other\") %>%",
    "    count(biological_program, recovery_direction, name = \"n_pathways\")",
    "",
    "  ggplot(",
    "    recovery_plot,",
    "    aes(",
    "      x = n_pathways,",
    "      y = fct_reorder(biological_program, n_pathways, .fun = sum),",
    "      fill = recovery_direction",
    "    )",
    "  ) +",
    "    geom_col() +",
    "    labs(",
    "      title = \"Recovery-direction pathways removed in Step 1\",",
    "      subtitle = \"Treatment FDR significant; Disease and Treatment NES have opposite signs.\",",
    "      x = \"Pathway rows\",",
    "      y = NULL,",
    "      fill = \"Recovery direction\"",
    "    ) +",
    "    theme_bw(base_size = 12) +",
    "    theme(",
    "      legend.position = \"bottom\",",
    "      axis.text.y = element_text(size = 10)",
    "    )",
    "}",
    "```",
    "",
    "```{r, results='asis'}",
    "if (nrow(recovery) > 0) {",
    "  recovery_summary <- recovery %>%",
    "    filter(biological_program != \"Other\") %>%",
    "    count(biological_program, sort = TRUE)",
    "",
    "  top_recovery <- recovery_summary %>%",
    "    slice_head(n = min(3, nrow(recovery_summary))) %>%",
    "    transmute(txt = paste0(biological_program, \" (n=\", n, \")\")) %>%",
    "    pull(txt)",
    "",
    "  ifng_n <- recovery_summary %>%",
    "    filter(biological_program == \"IFN-gamma / type II interferon signaling\") %>%",
    "    pull(n)",
    "",
    "  extra_ifng <- if (length(ifng_n) == 1) {",
    "    paste0(",
    "      \" IFN-gamma/type II interferon contributes \", ifng_n,",
    "      \" recovery-direction rows, which is useful as a biological calibration that the recovery filter is capturing disease-high pathways that decrease after treatment.\"",
    "    )",
    "  } else {",
    "    \"\"",
    "  }",
    "",
    "  comment_box(",
    "    paste0(",
    "      \"The largest recovery-associated biological programs are \",",
    "      paste(top_recovery, collapse = \"; \"), \".\",",
    "      extra_ifng,",
    "      \" Counts are pathway rows and therefore reflect overlapping pathway definitions rather than independent biological effects.\"",
    "    )",
    "  )",
    "}",
    "```",
    "",
    "## Step 1 / Step 2 / Step 3 contribution blocks",
    "",
    "This view shows how the **treatment-significant starting universe** is partitioned by the three-stage workflow. Because Step 1, Step 2, and Step 3 are nested, the plotted blocks are made mutually exclusive:",
    "",
    "- **Recovery removed before Step 1** — treatment-significant pathways with opposite Disease/Treatment NES directions.",
    "- **Step 1 only** — survives recovery removal but does not pass the Step 2 Disease-FDR criterion.",
    "- **Step 2 only** — passes Step 2 but does not pass the Step 3 NES prioritization.",
    "- **Step 3 retained** — survives all filters and represents the final prioritized set used for Questions 2 and 3.",
    "",
    "```{r}",
    "#| fig-width: 16",
    "#| fig-height: 10",
    "if (!is.null(step_contribution_blocks) && nrow(step_contribution_blocks) > 0) {",
    "  step_block_plot_df <- step_contribution_blocks %>%",
    "    count(database, celltype, step_block, name = \"n_pathways\") %>%",
    "    group_by(database, celltype) %>%",
    "    mutate(total = sum(n_pathways)) %>%",
    "    ungroup()",
    "",
    "  cell_order <- step_block_plot_df %>%",
    "    group_by(celltype) %>%",
    "    summarise(total_all = sum(n_pathways), .groups = \"drop\") %>%",
    "    arrange(total_all) %>%",
    "    pull(celltype)",
    "",
    "  step_block_plot_df <- step_block_plot_df %>%",
    "    mutate(celltype = factor(celltype, levels = cell_order))",
    "",
    "  ggplot(",
    "    step_block_plot_df,",
    "    aes(",
    "      x = n_pathways,",
    "      y = celltype,",
    "      fill = step_block",
    "    )",
    "  ) +",
    "    geom_col(width = 0.78) +",
    "    facet_wrap(~ database, scales = \"free_x\", nrow = 1) +",
    "    labs(",
    "      title = \"Contribution of Step 1, Step 2, and Step 3 to the treatment-significant pathway universe\",",
    "      subtitle = paste0(",
    "        \"Treatment FDR < \", settings$treat_fdr,",
    "        \"; Step 2: Disease FDR >= \", settings$disease_fdr_min,",
    "        \"; Step 3: |Treatment NES| >= \", settings$treatment_nes_min,",
    "        if (isTRUE(settings$use_disease_nes)) {",
    "          paste0(\" and |Disease NES| <= \", settings$disease_nes_max)",
    "        } else {",
    "          \"\"",
    "        },",
    "        \".\"",
    "      ),",
    "      x = \"Number of treatment-significant pathway rows\",",
    "      y = NULL,",
    "      fill = \"Analysis block\"",
    "    ) +",
    "    theme_bw(base_size = 12) +",
    "    theme(",
    "      legend.position = \"bottom\",",
    "      strip.text = element_text(face = \"bold\"),",
    "      axis.text.y = element_text(size = 9)",
    "    )",
    "}",
    "```",
    "",
    "```{r, results='asis'}",
    "if (!is.null(step_contribution_blocks) && nrow(step_contribution_blocks) > 0) {",
    "  block_summary <- step_contribution_blocks %>%",
    "    count(step_block, name = \"n_pathways\") %>%",
    "    mutate(percent = 100 * n_pathways / sum(n_pathways))",
    "",
    "  n_total <- sum(block_summary$n_pathways)",
    "",
    "  get_block_n <- function(label) {",
    "    z <- block_summary %>%",
    "      filter(as.character(step_block) == label) %>%",
    "      pull(n_pathways)",
    "    if (length(z) == 0) 0 else z[[1]]",
    "  }",
    "",
    "  n_recovery <- get_block_n(\"Recovery removed before Step 1\")",
    "  n_step1only <- get_block_n(\"Step 1 only\")",
    "  n_step2only <- get_block_n(\"Step 2 only\")",
    "  n_step3 <- get_block_n(\"Step 3 retained\")",
    "",
    "  comment_box(",
    "    paste0(",
    "      \"The treatment-significant starting universe contains \", fmt_n(n_total), \" pathway rows. \",",
    "      fmt_n(n_recovery), \" are removed as recovery before Step 1; \",",
    "      fmt_n(n_step1only), \" survive Step 1 but do not enter Step 2; \",",
    "      fmt_n(n_step2only), \" pass Step 2 but do not meet Step 3 NES prioritization; and \",",
    "      fmt_n(n_step3), \" remain in the final Step 3 prioritized set. \",",
    "      \"This decomposition shows which cell types and databases contribute most strongly at each stage of the workflow.\"",
    "    )",
    "  )",
    "}",
    "```",
    "",
    "### Step contribution by cell type and database",
    "",
    "```{r}",
    "if (!is.null(step_contribution_blocks) && nrow(step_contribution_blocks) > 0) {",
    "  step_block_table <- step_contribution_blocks %>%",
    "    count(database, celltype, step_block, name = \"n_pathways\") %>%",
    "    group_by(database, celltype) %>%",
    "    mutate(",
    "      total_treatment_significant = sum(n_pathways),",
    "      percent_of_celltype_database =",
    "        100 * n_pathways / total_treatment_significant",
    "    ) %>%",
    "    ungroup() %>%",
    "    arrange(database, desc(total_treatment_significant), celltype, step_block)",
    "",
    "  report_dt(",
    "    step_block_table,",
    "    page_length = 25,",
    "    filter_top = TRUE,",
    "    round_digits = 2",
    "  )",
    "}",
    "```",
    "",
    "### Overall contribution of each step block",
    "",
    "```{r}",
    "if (!is.null(step_contribution_blocks) && nrow(step_contribution_blocks) > 0) {",
    "  step_block_summary <- step_contribution_blocks %>%",
    "    count(step_block, name = \"n_pathways\") %>%",
    "    mutate(percent = 100 * n_pathways / sum(n_pathways))",
    "",
    "  report_dt(",
    "    step_block_summary,",
    "    page_length = 10,",
    "    filter_top = FALSE,",
    "    round_digits = 2",
    "  )",
    "}",
    "```",
    "",
    "## Step 1 survivors",
    "",
    "```{r}",
    "#| fig-width: 13",
    "#| fig-height: 7",
    "step1_counts <- step1 %>%",
    "  count(celltype, treatment_direction, name = \"n_pathways\")",
    "",
    "if (nrow(step1_counts) > 0) {",
    "  ggplot(",
    "    step1_counts,",
    "    aes(",
    "      x = n_pathways,",
    "      y = fct_reorder(celltype, n_pathways, .fun = sum),",
    "      fill = treatment_direction",
    "    )",
    "  ) +",
    "    geom_col() +",
    "    labs(",
    "      title = \"Step 1: treatment-responsive pathways after recovery removal\",",
    "      x = \"Pathway rows\",",
    "      y = NULL,",
    "      fill = NULL",
    "    ) +",
    "    theme_bw(base_size = 12) +",
    "    theme(legend.position = \"bottom\")",
    "}",
    "```",
    "",
    "# Step 2 — FDR-defined candidate set",
    "",
    "> **Question 1:** Which PBMC pathways change after treatment but are not significantly changed in untreated APECED versus healthy controls?",
    "",
    "Starting from Step 1 survivors, retain pathways with:",
    "",
    "- **Treatment FDR < selected cutoff**",
    "- **Disease FDR >= selected cutoff**",
    "",
    "NES magnitude is not used to determine membership in Step 2. This defines the formal treatment-associated, non-recovery candidate set used to answer Question 1. Disease FDR above the cutoff is treated as a lack of statistical significance in this small cohort, not as evidence that the disease effect is absent.",
    "",
    "```{r}",
    "#| fig-width: 13",
    "#| fig-height: 7",
    "step2_counts <- step2 %>%",
    "  count(celltype, treatment_direction, name = \"n_pathways\")",
    "",
    "if (nrow(step2_counts) > 0) {",
    "  ggplot(",
    "    step2_counts,",
    "    aes(",
    "      x = n_pathways,",
    "      y = fct_reorder(celltype, n_pathways, .fun = sum),",
    "      fill = treatment_direction",
    "    )",
    "  ) +",
    "    geom_col() +",
    "    labs(",
    "      title = \"Step 2: FDR-defined treatment-not-disease candidate set\",",
    "      x = \"Pathway rows\",",
    "      y = NULL,",
    "      fill = NULL",
    "    ) +",
    "    theme_bw(base_size = 12) +",
    "    theme(legend.position = \"bottom\")",
    "}",
    "```",
    "",
    "## Step 2 biological programs",
    "",
    "```{r}",
    "step2_programs <- step2 %>%",
    "  filter(biological_program != \"Other\") %>%",
    "  group_by(celltype, biological_program) %>%",
    "  summarise(",
    "    n_pathways = n(),",
    "    n_databases = n_distinct(database),",
    "    median_treatment_NES = median(treatment_NES, na.rm = TRUE),",
    "    median_abs_treatment_NES = median(abs_treatment_NES, na.rm = TRUE),",
    "    median_disease_NES = median(disease_NES, na.rm = TRUE),",
    "    best_treatment_FDR = min(treatment_FDR, na.rm = TRUE),",
    "    .groups = \"drop\"",
    "  ) %>%",
    "  arrange(desc(n_pathways), desc(n_databases), best_treatment_FDR)",
    "```",
    "",
    "```{r}",
    "#| fig-width: 17",
    "#| fig-height: 8",
    "if (nrow(step2_programs) > 0) {",
    "  ggplot(",
    "    step2_programs,",
    "    aes(",
    "      x = celltype,",
    "      y = biological_program,",
    "      size = n_pathways,",
    "      shape = n_databases >= 2",
    "    )",
    "  ) +",
    "    geom_point(alpha = 0.85, stroke = 0.8) +",
    "    scale_shape_manual(values = c(`FALSE` = 21, `TRUE` = 24)) +",
    "    labs(",
    "      title = \"Step 2 biological programs by PBMC cell type\",",
    "      subtitle = \"Size = retained pathway rows; triangle = biological program represented in at least two pathway databases.\",",
    "      x = \"Cell type\",",
    "      y = NULL,",
    "      size = \"Pathway rows\",",
    "      shape = \">=2 databases\"",
    "    ) +",
    "    theme_bw(base_size = 12) +",
    "    theme(",
    "      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 9),",
    "      axis.text.y = element_text(size = 10),",
    "      legend.position = \"bottom\",",
    "      plot.margin = margin(8, 15, 8, 8)",
    "    )",
    "}",
    "```",
    "",
    "```{r, results='asis'}",
    "if (nrow(step2_programs) > 0) {",
    "  step2_overall <- step2 %>%",
    "    filter(biological_program != \"Other\") %>%",
    "    group_by(biological_program) %>%",
    "    summarise(",
    "      n_pathways = n(),",
    "      n_celltypes = n_distinct(celltype),",
    "      n_databases = n_distinct(database),",
    "      median_abs_treatment_NES = median(abs_treatment_NES, na.rm = TRUE),",
    "      .groups = \"drop\"",
    "    ) %>%",
    "    arrange(desc(n_pathways), desc(n_celltypes), desc(n_databases))",
    "",
    "  top2 <- step2_overall %>%",
    "    slice_head(n = min(4, nrow(step2_overall))) %>%",
    "    transmute(",
    "      txt = paste0(",
    "        biological_program, \" (\", n_pathways, \" rows; \",",
    "        n_celltypes, \" cell types; \", n_databases, \" database\",",
    "        ifelse(n_databases == 1, \"\", \"s\"), \")\"",
    "      )",
    "    ) %>%",
    "    pull(txt)",
    "",
    "  comment_box(",
    "    paste0(",
    "      \"After recovery removal and the disease-FDR screen, the most represented Step 2 programs are \",",
    "      paste(top2, collapse = \"; \"), \". \",",
    "      \"Programs appearing across several cell types are more compelling as broad treatment-associated effects, and support from more than one pathway database strengthens biological-program consistency. Raw pathway-row count should be treated as secondary because program categories differ in the number of overlapping pathway definitions they contain.\"",
    "    )",
    "  )",
    "}",
    "```",
    "",
    "### Step 2 program table",
    "",
    "```{r}",
    "report_dt(step2_programs, page_length = 25, filter_top = TRUE, round_digits = 3)",
    "```",
    "",
    "",
    "## Answer to Question 1",
    "",
    "```{r, results='asis'}",
    "if (nrow(step2) == 0) {",
    "  comment_box(",
    "    \"No pathways meet the current Question 1 definition after recovery removal and the selected FDR thresholds.\"",
    "  )",
    "} else {",
    "  q1_programs <- step2 %>%",
    "    filter(biological_program != \"Other\") %>%",
    "    group_by(biological_program) %>%",
    "    summarise(",
    "      n_pathways = n(),",
    "      n_celltypes = n_distinct(celltype),",
    "      n_databases = n_distinct(database),",
    "      .groups = \"drop\"",
    "    ) %>%",
    "    arrange(desc(n_pathways), desc(n_celltypes), desc(n_databases))",
    "",
    "  q1_top <- q1_programs %>%",
    "    slice_head(n = min(5, nrow(q1_programs))) %>%",
    "    transmute(",
    "      txt = paste0(",
    "        biological_program,",
    "        \" (\", n_pathways, \" rows; \",",
    "        n_celltypes, \" cell types; \",",
    "        n_databases, \" database\",",
    "        ifelse(n_databases == 1, \"\", \"s\"), \")\"",
    "      )",
    "    ) %>%",
    "    pull(txt)",
    "",
    "  comment_box(",
    "    paste0(",
    "      \"Question 1 identifies \", fmt_n(nrow(step2)),",
    "      \" pathway rows that are significantly treatment-responsive, survive recovery removal, and are not statistically significant in the untreated-versus-healthy PBMC comparison. \",",
    "      \"The most represented biological programs are \",",
    "      paste(q1_top, collapse = \"; \"), \". \",",
    "      \"These are treatment-associated, non-recovery candidates; disease FDR above the cutoff should be interpreted as non-significant in this dataset, not as proof of no disease effect.\"",
    "    )",
    "  )",
    "}",
    "```",
    "",
    "# Step 3 — NES prioritization",
    "",
    "> **Question 2:** Which of the Question 1 candidates show the strongest treatment effects and the strongest biological-program support across pathway databases?",
    "",
    "Step 3 prioritizes the strongest Step 2 candidates using treatment-effect magnitude and, when enabled, a maximum residual disease NES.",
    "",
    "The thresholds shown below are the exact values selected in the app for this report. The optional Disease NES threshold is used as an effect-size prioritization tool rather than an equivalence test.",
    "",
    "```{r}",
    "step3_programs <- step3 %>%",
    "  filter(biological_program != \"Other\") %>%",
    "  group_by(celltype, biological_program) %>%",
    "  summarise(",
    "    n_pathways = n(),",
    "    n_databases = n_distinct(database),",
    "    n_GO_BP = sum(database == \"GO_BP\"),",
    "    n_Reactome = sum(database == \"Reactome\"),",
    "    n_Hallmark = sum(database == \"Hallmark\"),",
    "    median_treatment_NES = median(treatment_NES, na.rm = TRUE),",
    "    median_abs_treatment_NES = median(abs_treatment_NES, na.rm = TRUE),",
    "    median_disease_NES = median(disease_NES, na.rm = TRUE),",
    "    median_abs_disease_NES = median(abs_disease_NES, na.rm = TRUE),",
    "    best_treatment_FDR = min(treatment_FDR, na.rm = TRUE),",
    "    .groups = \"drop\"",
    "  ) %>%",
    "  arrange(",
    "    desc(n_databases),",
    "    desc(n_pathways),",
    "    best_treatment_FDR,",
    "    desc(median_abs_treatment_NES)",
    "  )",
    "```",
    "",
    "## Step 3 biological programs",
    "",
    "```{r}",
    "#| fig-width: 17",
    "#| fig-height: 8.5",
    "if (nrow(step3_programs) > 0) {",
    "",
    "  # Make the biologically important feature — breadth across PBMC cell types —",
    "  # explicit without changing which pathways pass the filters.",
    "  program_breadth <- step3_programs %>%",
    "    group_by(biological_program) %>%",
    "    summarise(",
    "      n_celltypes = n_distinct(celltype),",
    "      total_pathway_rows = sum(n_pathways),",
    "      overall_median_treatment_NES = weighted.mean(",
    "        median_treatment_NES,",
    "        w = pmax(n_pathways, 1),",
    "        na.rm = TRUE",
    "      ),",
    "      .groups = \"drop\"",
    "    ) %>%",
    "    arrange(desc(n_celltypes), desc(abs(overall_median_treatment_NES)), desc(total_pathway_rows))",
    "",
    "  total_selected_celltypes <- length(unique(step3$celltype))",
    "",
    "  program_breadth <- program_breadth %>%",
    "    mutate(",
    "      breadth_label = paste0(",
    "        biological_program,",
    "        \"  [\", n_celltypes, \"/\", total_selected_celltypes, \" cell types]\"",
    "      )",
    "    )",
    "",
    "  plot_df <- step3_programs %>%",
    "    left_join(",
    "      program_breadth %>% select(biological_program, n_celltypes, breadth_label),",
    "      by = \"biological_program\"",
    "    ) %>%",
    "    mutate(",
    "      breadth_label = factor(",
    "        breadth_label,",
    "        levels = rev(program_breadth$breadth_label)",
    "      )",
    "    )",
    "",
    "  # Keep a symmetric NES color scale centered at zero so up/down effects",
    "  # are visually comparable.",
    "  nes_limit <- max(abs(plot_df$median_treatment_NES), na.rm = TRUE)",
    "  if (!is.finite(nes_limit) || nes_limit == 0) nes_limit <- 1",
    "",
    "  ggplot(",
    "    plot_df,",
    "    aes(",
    "      x = celltype,",
    "      y = breadth_label,",
    "      size = n_pathways,",
    "      fill = median_treatment_NES,",
    "      shape = n_databases >= 2",
    "    )",
    "  ) +",
    "    geom_point(alpha = 0.95, stroke = 0.9, colour = \"grey20\") +",
    "    scale_size_area(max_size = 11) +",
    "    scale_fill_gradient2(",
    "      low = \"#2166AC\",",
    "      mid = \"white\",",
    "      high = \"#B2182B\",",
    "      midpoint = 0,",
    "      limits = c(-nes_limit, nes_limit),",
    "      oob = scales::squish",
    "    ) +",
    "    scale_shape_manual(values = c(`FALSE` = 21, `TRUE` = 24)) +",
    "    labs(",
    "      title = \"Step 3: strongest remaining biological programs by PBMC cell type\",",
    "      subtitle = paste0(",
    "        \"Programs are ordered by PBMC cell-type breadth. \",",
    "        \"Y-axis labels show the number of cell types with a retained signal. \",",
    "        \"Size = pathway rows; fill = median Treatment NES; triangle = support from >=2 databases.\"",
    "      ),",
    "      x = \"Cell type\",",
    "      y = NULL,",
    "      size = \"Pathway rows\",",
    "      fill = \"Median Treatment NES\",",
    "      shape = \">=2 databases\"",
    "    ) +",
    "    theme_bw(base_size = 12) +",
    "    theme(",
    "      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 9),",
    "      axis.text.y = element_text(size = 10.5, face = \"bold\"),",
    "      legend.position = \"bottom\",",
    "      panel.grid.minor = element_blank(),",
    "      plot.margin = margin(8, 20, 8, 8)",
    "    )",
    "}",
    "```",
    "",
    "```{r, results='asis'}",
    "if (nrow(step3_programs) > 0) {",
    "  step3_overall <- step3 %>%",
    "    filter(biological_program != \"Other\") %>%",
    "    group_by(biological_program) %>%",
    "    summarise(",
    "      n_pathways = n(),",
    "      n_celltypes = n_distinct(celltype),",
    "      n_databases = n_distinct(database),",
    "      median_treatment_NES = median(treatment_NES, na.rm = TRUE),",
    "      median_abs_treatment_NES = median(abs(treatment_NES), na.rm = TRUE),",
    "      median_abs_disease_NES = median(abs(disease_NES), na.rm = TRUE),",
    "      best_treatment_FDR = min(treatment_FDR, na.rm = TRUE),",
    "      .groups = \"drop\"",
    "    ) %>%",
    "    arrange(",
    "      desc(n_databases),",
    "      desc(n_celltypes),",
    "      desc(n_pathways),",
    "      desc(median_abs_treatment_NES)",
    "    )",
    "",
    "  top3 <- step3_overall %>%",
    "    slice_head(n = min(5, nrow(step3_overall))) %>%",
    "    transmute(",
    "      txt = paste0(",
    "        biological_program,",
    "        \" (\", n_pathways, \" rows; \",",
    "        n_celltypes, \" cell types; \",",
    "        n_databases, \" database\", ifelse(n_databases == 1, \"\", \"s\"),",
    "        \"; median |Treatment NES|=\", sprintf(\"%.2f\", median_abs_treatment_NES), \")\"",
    "      )",
    "    ) %>%",
    "    pull(txt)",
    "",
    "  direction_summary <- step3 %>%",
    "    filter(biological_program != \"Other\") %>%",
    "    group_by(biological_program) %>%",
    "    summarise(",
    "      n_pathways = n(),",
    "      median_treatment_NES = median(treatment_NES, na.rm = TRUE),",
    "      n_positive = sum(treatment_NES > 0, na.rm = TRUE),",
    "      n_negative = sum(treatment_NES < 0, na.rm = TRUE),",
    "      .groups = \"drop\"",
    "    ) %>%",
    "    mutate(",
    "      dominant_direction = case_when(",
    "        n_positive > n_negative ~ \"predominantly increased\",",
    "        n_negative > n_positive ~ \"predominantly decreased\",",
    "        TRUE ~ \"mixed direction\"",
    "      )",
    "    ) %>%",
    "    arrange(desc(n_pathways))",
    "",
    "  dir_top <- direction_summary %>%",
    "    slice_head(n = min(5, nrow(direction_summary))) %>%",
    "    transmute(",
    "      txt = paste0(",
    "        biological_program,",
    "        \" (\", n_pathways, \" rows; \",",
    "        dominant_direction,",
    "        \"; median Treatment NES=\",",
    "        sprintf(\"%.2f\", median_treatment_NES), \")\"",
    "      )",
    "    ) %>%",
    "    pull(txt)",
    "",
    "  comment_box(",
    "    paste0(",
    "      \"Question 2 prioritizes the strongest remaining treatment-associated programs primarily by breadth across PBMC cell types, then by Treatment NES direction/magnitude and cross-database support. Pathway-set size is treated as secondary because program categories contain unequal numbers of overlapping pathway definitions. \",",
    "      \"The leading programs are \",",
    "      paste(dir_top, collapse = \"; \"), \". \",",
    "      \"Larger bubbles indicate more retained pathway rows, while the median Treatment NES indicates whether that biological program is generally increased or decreased after treatment. \",",
    "      \"Support from at least two pathway databases strengthens biological-program consistency but is not exact pathway replication.\"",
    "    )",
    "  )",
    "}",
    "```",
    "",
    "",
    "## Interpretation of the Step 3 pattern",
    "",
    "The Step 3 plot is also the key visual summary for identifying **candidate broader pharmacologic effects** that remain after disease-recovery pathways have been removed.",
    "",
    "For biological interpretation, **cell-type breadth is treated as a primary metric**: a program observed across many distinct PBMC cell types is more consistent with a broad systemic treatment-associated effect than one confined to a small number of lineages. Raw pathway-row count is retained as a secondary measure of within-program support because pathway databases differ in how many overlapping pathway definitions they contain.",
    "",
    "```{r, results='asis'}",
    "if (nrow(step3) > 0) {",
    "  residual_programs <- step3 %>%",
    "    filter(biological_program != \"Other\") %>%",
    "    group_by(biological_program) %>%",
    "    summarise(",
    "      n_pathways = n(),",
    "      n_celltypes = n_distinct(celltype),",
    "      n_databases = n_distinct(database),",
    "      median_treatment_NES = median(treatment_NES, na.rm = TRUE),",
    "      .groups = \"drop\"",
    "    ) %>%",
    "    arrange(desc(n_pathways), desc(n_celltypes))",
    "",
    "  get_prog <- function(label) {",
    "    residual_programs %>% filter(biological_program == label)",
    "  }",
    "",
    "  tr <- get_prog(\"Translation / ribosome / RNA quality control\")",
    "  cc <- get_prog(\"Cell cycle / proliferation\")",
    "  ii <- get_prog(\"Other immune / innate signaling\")",
    "  js <- get_prog(\"Other cytokine / JAK-STAT signaling\")",
    "  ig <- get_prog(\"IFN-gamma / type II interferon signaling\")",
    "  ia <- get_prog(\"IFN-alpha/beta / type I interferon signaling\")",
    "",
    "  recovery_ifng_n <- recovery %>%",
    "    filter(biological_program == \"IFN-gamma / type II interferon signaling\") %>%",
    "    nrow()",
    "",
    "  residual_ifng_n <- if (nrow(ig) == 1) ig$n_pathways[[1]] else 0",
    "",
    "  parts <- character()",
    "",
    "  if (nrow(tr) == 1) {",
    "    parts <- c(parts, paste0(",
    "      \"Translation/ribosome/RNA-quality-control is a broad residual treatment-associated program (\",",
    "      tr$n_pathways[[1]], \" pathway rows across \", tr$n_celltypes[[1]],",
    "      \" PBMC cell types; median Treatment NES=\", sprintf(\"%.2f\", tr$median_treatment_NES[[1]]), \" ).\"",
    "    ))",
    "  }",
    "",
    "  if (nrow(cc) == 1) {",
    "    parts <- c(parts, paste0(",
    "      \"Cell-cycle/proliferation is also prominent (\",",
    "      cc$n_pathways[[1]], \" pathway rows across \", cc$n_celltypes[[1]],",
    "      \" cell types; median Treatment NES=\", sprintf(\"%.2f\", cc$median_treatment_NES[[1]]), \" ).\"",
    "    ))",
    "  }",
    "",
    "  if (nrow(ii) == 1 && is.finite(ii$median_treatment_NES[[1]]) && ii$median_treatment_NES[[1]] < 0) {",
    "    parts <- c(parts, paste0(",
    "      \"Residual immune/innate signaling is predominantly decreased after treatment (\",",
    "      ii$n_pathways[[1]], \" pathway rows; median Treatment NES=\",",
    "      sprintf(\"%.2f\", ii$median_treatment_NES[[1]]), \" ), consistent with broader immune suppression beyond pathways classified as disease recovery.\"",
    "    ))",
    "  }",
    "",
    "  if (nrow(js) == 1) {",
    "    parts <- c(parts, paste0(",
    "      \"Other cytokine/JAK-STAT signaling remains represented in the residual set, supporting the interpretation that some retained changes may reflect broader pharmacologic pathway inhibition rather than correction of APECED-associated biology.\"",
    "    ))",
    "  }",
    "",
    "  if (recovery_ifng_n > residual_ifng_n) {",
    "    parts <- c(parts, paste0(",
    "      \"IFN-gamma/type II interferon contributes more strongly to the recovery set than to the final residual set (\",",
    "      recovery_ifng_n, \" recovery-direction rows versus \", residual_ifng_n,",
    "      \" Step 3 rows). This provides a useful biological calibration that the recovery filter is removing expected disease-correction biology before broader treatment effects are interpreted.\"",
    "    ))",
    "  } else if (residual_ifng_n > 0) {",
    "    parts <- c(parts, paste0(",
    "      \"A smaller residual IFN-gamma/type II interferon component remains in Step 3 (\",",
    "      residual_ifng_n, \" pathway rows), so interferon biology is not completely absent from the residual treatment signature.\"",
    "    ))",
    "  }",
    "",
    "  type1_n <- if (nrow(ia) == 1) ia$n_pathways[[1]] else 0",
    "  if (type1_n > 0) {",
    "    parts <- c(parts, paste0(",
    "      \"Type I interferon signaling contributes \", type1_n,",
    "      \" residual pathway rows and is not the dominant Step 3 program.\"",
    "    ))",
    "  }",
    "",
    "  if (length(parts) > 0) {",
    "    comment_box(paste(parts, collapse = \" \"))",
    "  }",
    "}",
    "```",
    "",
    "> **Working conclusion:** After removal of pathways consistent with disease recovery, the strongest remaining treatment-associated effects are best interpreted as **candidate broader pharmacologic consequences of JAK/STAT inhibition**. Priority is given to programs that are **widespread across PBMC cell types**, show a substantial Treatment NES, and, where available, have support across pathway databases. Pathway-row count is used only as a secondary indicator because biological-program categories differ in the number of overlapping pathway definitions they contain. These findings should not be described as established adverse effects without independent functional or clinical validation.",
    "",
    "## Answer to Question 2",
    "",
    "The Step 3 biological-program plot is the primary summary for Question 2. The most important metric is **how broadly a program is distributed across PBMC cell types**. The **median Treatment NES** gives the direction and magnitude of the treatment effect, and representation in at least two databases provides stronger cross-database biological-program support. Bubble size (retained pathway rows) is informative but should be interpreted as a secondary measure because pathway databases contain unequal numbers of overlapping pathway definitions for different biological programs.",
    "",
    "## Strongest Step 3 programs",
    "",
    "```{r}",
    "if (nrow(step3_programs) == 0) {",
    "  cat(\"No biological programs pass the current Step 3 settings.\")",
    "} else {",
    "  report_dt(step3_programs, page_length = 25, filter_top = TRUE, round_digits = 3)",
    "}",
    "```",
    "",
    "## Strongest individual Step 3 pathways",
    "",
    "The table below adds **NES extremeness percentiles** for both comparisons. These are descriptive rankings, not additional significance tests. For each pathway, |NES| is ranked against **all pathways tested before filtering** within the same cell type and pathway database. Thus a Treatment |NES| percentile of 99 means that the treatment effect is more extreme than approximately 99% of pathways tested in that cell-type/database stratum.",
    "",
    "```{r}",
    "step3_top <- step3 %>%",
    "  arrange(",
    "    desc(treatment_abs_NES_percentile),",
    "    disease_abs_NES_percentile,",
    "    treatment_FDR,",
    "    desc(abs_treatment_NES)",
    "  ) %>%",
    "  select(",
    "    celltype,",
    "    database,",
    "    clean_pathway,",
    "    biological_program,",
    "    treatment_NES,",
    "    treatment_abs_NES_percentile,",
    "    treatment_FDR,",
    "    disease_NES,",
    "    disease_abs_NES_percentile,",
    "    disease_FDR,",
    "    extremeness_gap_percentile",
    "  ) %>%",
    "  rename(",
    "    `Treatment |NES| percentile` = treatment_abs_NES_percentile,",
    "    `Disease |NES| percentile` = disease_abs_NES_percentile,",
    "    `Treatment - disease percentile gap` = extremeness_gap_percentile",
    "  )",
    "",
    "if (nrow(step3_top) == 0) {",
    "  cat(\"No pathways pass the current Step 3 settings.\")",
    "} else {",
    "  report_dt(step3_top, page_length = 25, filter_top = TRUE, round_digits = 2)",
    "}",
    "```",
    "",
    "## Individual pathway dumbbell contrasts",
    "",
    "The dumbbell plots now follow the **top three biological programs from the Step 3 biological-program ranking shown above**. This keeps the individual-pathway figures aligned with the main biological conclusion rather than allowing a weaker program to dominate because of a few extreme pathway rows.",
    "",
    "For each of the top three Step 3 programs, the report creates a separate NES dumbbell plot showing the strongest individual pathways within that program. Pathways are ranked within each program by `|Treatment NES| - |Disease NES|`, with treatment FDR and |Treatment NES| used as tie-breakers. Point color identifies the comparison, while point shape shows FDR status: filled = FDR below the selected cutoff; open = FDR at or above the cutoff.",
    "",
    "```{r}",
    "# Use the same ranking logic as the Step 3 biological-program plot:",
    "# cell-type breadth first, then absolute overall treatment effect, then row count.",
    "dumbbell_program_rank <- step3 %>%",
    "  filter(",
    "    biological_program != \"Other\",",
    "    !is.na(treatment_NES),",
    "    !is.na(disease_NES),",
    "    !is.na(treatment_FDR),",
    "    !is.na(disease_FDR)",
    "  ) %>%",
    "  group_by(biological_program) %>%",
    "  summarise(",
    "    n_celltypes = n_distinct(celltype),",
    "    n_databases = n_distinct(database),",
    "    n_pathways = n(),",
    "    overall_median_treatment_NES = median(treatment_NES, na.rm = TRUE),",
    "    median_abs_treatment_NES = median(abs(treatment_NES), na.rm = TRUE),",
    "    .groups = \"drop\"",
    "  ) %>%",
    "  arrange(",
    "    desc(n_celltypes),",
    "    desc(abs(overall_median_treatment_NES)),",
    "    desc(n_pathways)",
    "  )",
    "",
    "top3_dumbbell_programs <- dumbbell_program_rank %>%",
    "  slice_head(n = 3) %>%",
    "  pull(biological_program)",
    "",
    "dumbbell_top3_df <- step3 %>%",
    "  filter(",
    "    biological_program %in% top3_dumbbell_programs,",
    "    !is.na(treatment_NES),",
    "    !is.na(disease_NES),",
    "    !is.na(treatment_FDR),",
    "    !is.na(disease_FDR)",
    "  ) %>%",
    "  mutate(",
    "    NES_abs_gap = abs(treatment_NES) - abs(disease_NES),",
    "    program_rank = match(biological_program, top3_dumbbell_programs)",
    "  ) %>%",
    "  group_by(biological_program) %>%",
    "  arrange(desc(NES_abs_gap), treatment_FDR, desc(abs(treatment_NES)), .by_group = TRUE) %>%",
    "  slice_head(n = 10) %>%",
    "  mutate(pathway_rank_within_program = row_number()) %>%",
    "  ungroup() %>%",
    "  arrange(program_rank, pathway_rank_within_program) %>%",
    "  mutate(",
    "    pathway_label = paste0(clean_pathway, \"  [\", celltype, \" | \", database, \"]\")",
    "  )",
    "```",
    "",
    "### Top three biological programs used for the dumbbell plots",
    "",
    "```{r}",
    "if (length(top3_dumbbell_programs) > 0) {",
    "  dumbbell_program_table <- dumbbell_program_rank %>%",
    "    filter(biological_program %in% top3_dumbbell_programs) %>%",
    "    mutate(program_rank = match(biological_program, top3_dumbbell_programs)) %>%",
    "    arrange(program_rank) %>%",
    "    select(program_rank, biological_program, n_celltypes, n_databases, n_pathways, overall_median_treatment_NES, median_abs_treatment_NES)",
    "",
    "  report_dt(dumbbell_program_table, page_length = 3, filter_top = FALSE, round_digits = 3)",
    "}",
    "```",
    "",
    "### Individual pathway NES dumbbells for the top three programs",
    "",
    "```{r}",
    "#| fig-width: 17",
    "#| fig-height: 8",
    "if (!requireNamespace(\"dumbbell\", quietly = TRUE)) {",
    "  cat(\"The R package 'dumbbell' is required for these figures. Install it with install.packages('dumbbell').\")",
    "} else if (length(top3_dumbbell_programs) == 0 || nrow(dumbbell_top3_df) == 0) {",
    "  cat(\"No Step 3 pathways with complete NES/FDR values are available for the top-three-program dumbbells.\")",
    "} else {",
    "  for (prog in top3_dumbbell_programs) {",
    "    prog_df <- dumbbell_top3_df %>%",
    "      filter(biological_program == prog) %>%",
    "      arrange(pathway_rank_within_program) %>%",
    "      mutate(",
    "        pathway_label = factor(pathway_label, levels = rev(unique(pathway_label))),",
    "        dumbbell_key = prog,",
    "        treatment_FDR_status = if_else(",
    "          treatment_FDR < settings$treat_fdr,",
    "          \"FDR < cutoff\",",
    "          \"FDR >= cutoff\"",
    "        ),",
    "        disease_FDR_status = if_else(",
    "          disease_FDR < settings$treat_fdr,",
    "          \"FDR < cutoff\",",
    "          \"FDR >= cutoff\"",
    "        )",
    "      )",
    "",
    "    if (nrow(prog_df) == 0) next",
    "",
    "    cat(\"\\n\\n#### \", prog, \"\\n\\n\", sep = \"\")",
    "",
    "    # Use dumbbell for the connecting segments, then draw the points separately",
    "    # so FDR significance can be encoded by point shape.",
    "    p_prog <- dumbbell::dumbbell(",
    "      xdf = prog_df,",
    "      id = \"pathway_label\",",
    "      key = \"dumbbell_key\",",
    "      column1 = \"disease_NES\",",
    "      column2 = \"treatment_NES\",",
    "      lab1 = \"Untreated vs Healthy\",",
    "      lab2 = \"Treated vs Untreated\",",
    "      pt_val = 0,",
    "      delt = 0,",
    "      textsize = 3,",
    "      pointsize = 0,",
    "      segsize = 1,",
    "      title = paste0(\"Individual pathway effect-size contrast — \", prog)",
    "    ) +",
    "      geom_point(",
    "        data = prog_df,",
    "        aes(x = disease_NES, y = pathway_label, shape = disease_FDR_status),",
    "        colour = \"blue\",",
    "        size = 3.4,",
    "        stroke = 1.1,",
    "        inherit.aes = FALSE",
    "      ) +",
    "      geom_point(",
    "        data = prog_df,",
    "        aes(x = treatment_NES, y = pathway_label, shape = treatment_FDR_status),",
    "        colour = \"red\",",
    "        size = 3.4,",
    "        stroke = 1.1,",
    "        inherit.aes = FALSE",
    "      ) +",
    "      geom_text(",
    "        data = prog_df,",
    "        aes(x = disease_NES, y = pathway_label, label = sprintf(\"%.2f\", disease_NES)),",
    "        nudge_y = -0.18,",
    "        size = 3,",
    "        colour = \"blue\",",
    "        inherit.aes = FALSE",
    "      ) +",
    "      geom_text(",
    "        data = prog_df,",
    "        aes(x = treatment_NES, y = pathway_label, label = sprintf(\"%.2f\", treatment_NES)),",
    "        nudge_y = -0.18,",
    "        size = 3,",
    "        colour = \"red\",",
    "        inherit.aes = FALSE",
    "      ) +",
    "      scale_shape_manual(",
    "        values = c(\"FDR < cutoff\" = 16, \"FDR >= cutoff\" = 1),",
    "        name = \"FDR significance\"",
    "      ) +",
    "      labs(",
    "        subtitle = paste0(",
    "          \"Top 10 pathways within this Step 3 biological program ranked by |Treatment NES| - |Disease NES|. \",",
    "          \"Red = Treated vs Untreated; blue = Untreated vs Healthy; filled = FDR < \",",
    "          settings$treat_fdr, \"; open = FDR >= \", settings$treat_fdr, \".\"",
    "        ),",
    "        x = \"Normalized enrichment score (NES)\",",
    "        y = NULL",
    "      ) +",
    "      geom_vline(xintercept = 0, linetype = \"dashed\", colour = \"grey55\") +",
    "      theme_bw(base_size = 11) +",
    "      theme(legend.position = \"bottom\", axis.text.y = element_text(size = 8.5))",
    "",
    "    print(p_prog)",
    "  }",
    "}",
    "```",
    "",
    "```{r, results='asis'}",
    "if (length(top3_dumbbell_programs) > 0) {",
    "  represented_programs <- dumbbell_program_rank %>%",
    "    filter(biological_program %in% top3_dumbbell_programs) %>%",
    "    mutate(program_rank = match(biological_program, top3_dumbbell_programs)) %>%",
    "    arrange(program_rank) %>%",
    "    transmute(",
    "      txt = paste0(",
    "        biological_program, \" (\", n_celltypes, \" cell types; \", n_databases, \" database\",",
    "        ifelse(n_databases == 1, \"\", \"s\"), \"; median |Treatment NES|=\",",
    "        sprintf(\"%.2f\", median_abs_treatment_NES), \")\"",
    "      )",
    "    ) %>%",
    "    pull(txt)",
    "",
    "  comment_box(",
    "    paste0(",
    "      \"The three dumbbell plots correspond directly to the top three biological programs in the Step 3 breadth ranking: \",",
    "      paste(represented_programs, collapse = \"; \"), \". \",",
    "      \"Each plot shows up to ten of the strongest individual pathways within that program, allowing pathway-level inspection without changing the program-level hierarchy used for Question 2.\"",
    "    )",
    "  )",
    "}",
    "```",
    "",
    "# Disease NES vs treatment NES after filtering",
    "",
    "```{r}",
    "#| fig-width: 16",
    "#| fig-height: 10",
    "if (nrow(step3) > 0) {",
    "  ggplot(",
    "    step3,",
    "    aes(",
    "      x = disease_NES,",
    "      y = treatment_NES,",
    "      shape = database",
    "    )",
    "  ) +",
    "    geom_hline(yintercept = 0, linetype = \"dashed\") +",
    "    geom_vline(xintercept = 0, linetype = \"dashed\") +",
    "    geom_point(alpha = 0.7, size = 2) +",
    "    facet_wrap(~celltype, scales = \"free\", ncol = 4) +",
    "    labs(",
    "      title = \"Step 3 disease NES vs treatment NES\",",
    "      subtitle = \"Only pathways surviving all three analysis steps are shown.\",",
    "      x = \"Untreated APECED vs healthy NES\",",
    "      y = \"Treated vs untreated NES\",",
    "      shape = \"Database\"",
    "    ) +",
    "    theme_bw(base_size = 11) +",
    "    theme(",
    "      legend.position = \"bottom\",",
    "      strip.text = element_text(face = \"bold\", size = 9)",
    "    )",
    "}",
    "```",
    "",
    "```{r, results='asis'}",
    "if (nrow(step3) > 0) {",
    "  n_same_pos <- sum(step3$disease_NES > 0 & step3$treatment_NES > 0, na.rm = TRUE)",
    "  n_same_neg <- sum(step3$disease_NES < 0 & step3$treatment_NES < 0, na.rm = TRUE)",
    "  n_missing <- sum(is.na(step3$disease_NES) | is.na(step3$treatment_NES))",
    "",
    "  comment_box(",
    "    paste0(",
    "      \"Because opposite-direction recovery rows were removed in Step 1, the remaining Step 3 points should predominantly represent same-direction or weak/near-zero disease effects rather than classical recovery. \",",
    "      \"Among rows with non-missing NES values, \", n_same_pos,",
    "      \" are disease-positive/treatment-positive and \", n_same_neg,",
    "      \" are disease-negative/treatment-negative. \",",
    "      if (n_missing > 0) paste0(n_missing, \" rows have missing NES values. \") else \"\",",
    "      \"Direction should be interpreted together with FDR and NES magnitude rather than as evidence of equivalence to healthy controls.\"",
    "    )",
    "  )",
    "}",
    "```",
    "",
    "# Question 3 — Candidate broader pharmacologic effects",
    "",
    "> **Question 3:** Which of the strongest non-recovery treatment-associated effects could represent broader pharmacologic consequences of JAK/STAT inhibition rather than correction of APECED-associated biology?",
    "",
    "The Step 3 results are the preferred set for this biological interpretation because these pathways have:",
    "",
    "1. survived removal of treatment responses consistent with disease recovery;",
    "2. shown a statistically significant treatment response while not reaching significance in the disease comparison; and",
    "3. passed the selected NES prioritization thresholds.",
    "",
    "```{r, results='asis'}",
    "if (nrow(step3) == 0) {",
    "  comment_box(",
    "    \"No pathways remain for Question 3 interpretation under the current Step 3 settings.\"",
    "  )",
    "} else {",
    "  q3_programs <- step3 %>%",
    "    filter(biological_program != \"Other\") %>%",
    "    group_by(biological_program) %>%",
    "    summarise(",
    "      n_pathways = n(),",
    "      n_celltypes = n_distinct(celltype),",
    "      n_databases = n_distinct(database),",
    "      median_treatment_NES = median(treatment_NES, na.rm = TRUE),",
    "      n_positive = sum(treatment_NES > 0, na.rm = TRUE),",
    "      n_negative = sum(treatment_NES < 0, na.rm = TRUE),",
    "      .groups = \"drop\"",
    "    ) %>%",
    "    mutate(",
    "      direction = case_when(",
    "        n_positive > n_negative ~ \"increased after treatment\",",
    "        n_negative > n_positive ~ \"decreased after treatment\",",
    "        TRUE ~ \"mixed direction\"",
    "      )",
    "    ) %>%",
    "    arrange(",
    "      desc(n_databases),",
    "      desc(n_celltypes),",
    "      desc(n_pathways)",
    "    )",
    "",
    "  q3_top <- q3_programs %>%",
    "    slice_head(n = min(6, nrow(q3_programs))) %>%",
    "    transmute(",
    "      txt = paste0(",
    "        biological_program,",
    "        \" (\", n_pathways, \" rows; \",",
    "        n_celltypes, \" cell types; \",",
    "        n_databases, \" database\",",
    "        ifelse(n_databases == 1, \"\", \"s\"),",
    "        \"; \", direction, \")\"",
    "      )",
    "    ) %>%",
    "    pull(txt)",
    "",
    "  comment_box(",
    "    paste0(",
    "      \"The leading Question 3 candidates are \",",
    "      paste(q3_top, collapse = \"; \"), \". \",",
    "      \"These are the strongest candidates for broader pharmacologic effects because they remain after recovery removal and Step 2/3 prioritization. \",",
    "      \"They should still be described as candidate broader treatment-associated effects rather than established adverse effects.\"",
    "    )",
    "  )",
    "}",
    "```",
    "",
    "### Overall interpretation",
    "",
    "The central result is that the residual Step 3 signal is not dominated by the interferon programs expected to reflect disease correction. Instead, several non-interferon programs remain distributed across multiple PBMC cell types, including translation/ribosome/RNA-quality-control, cell-cycle/proliferation, and residual immune/innate or cytokine/JAK-STAT effects. Their **cell-type breadth**, together with Treatment NES magnitude and cross-database support where available, is more informative than raw pathway-row count alone. This separation supports the working hypothesis that the analysis is distinguishing **disease recovery** from **broader treatment-associated biology**.",
    "",
    "These residual programs are therefore candidates for follow-up as potentially avoidable pharmacologic effects of broad JAK/STAT inhibition that might be reduced with a more targeted therapeutic strategy. This remains a hypothesis-generating interpretation rather than evidence of clinical toxicity.",
    "",
    "A disease FDR above the cutoff means **not statistically significant in this PBMC comparison**; it is not evidence that the disease effect is exactly zero. Likewise, GO BP, Reactome, and Hallmark contain overlapping pathway definitions, so repeated rows should be interpreted as biological-program consistency rather than independent replication.",
    "",
    "# Appendix — Step 1 recovery pathways",
    "",
    "```{r}",
    "if (nrow(recovery) > 0) {",
    "  recovery_appendix <- recovery %>%",
    "    select(",
    "      celltype,",
    "      database,",
    "      clean_pathway,",
    "      biological_program,",
    "      recovery_direction,",
    "      treatment_NES,",
    "      treatment_FDR,",
    "      disease_NES,",
    "      disease_FDR",
    "    )",
    "",
    "  report_dt(recovery_appendix, page_length = 25, filter_top = TRUE, round_digits = 4)",
    "}",
    "```",
    "",
    "# Appendix — Step 2 candidate pathways",
    "",
    "```{r}",
    "if (nrow(step2) > 0) {",
    "  step2_appendix <- step2 %>%",
    "    select(",
    "      celltype,",
    "      database,",
    "      clean_pathway,",
    "      biological_program,",
    "      treatment_direction,",
    "      treatment_NES,",
    "      treatment_FDR,",
    "      disease_NES,",
    "      disease_FDR",
    "    )",
    "",
    "  report_dt(step2_appendix, page_length = 25, filter_top = TRUE, round_digits = 4)",
    "}",
    "```"
  )
  }

  # -----------------------------------------------------------------------
  # REPORT SNAPSHOTS
  # -----------------------------------------------------------------------
  # These objects reproduce all three stages regardless of which step is
  # currently displayed in the app.
  report_step1 <- reactive({
    recovery_cleaned_treatment_universe()
  })

  report_step2 <- reactive({
    report_step1() %>%
      filter(
        !is.na(disease_FDR),
        disease_FDR >= input$disease_fdr_min
      )
  })


  # -----------------------------------------------------------------------
  # INTERFERON CLASSIFICATION QC
  # -----------------------------------------------------------------------
  # Shows every interferon-labeled pathway in the selected cell types/databases,
  # its manual biological-program assignment, and where it lands in the
  # recovery -> Step 1 -> Step 2 -> Step 3 workflow.
  interferon_classification_qc <- reactive({
    req(length(input$db) > 0)
    req(length(input$celltypes) > 0)

    ifn_base <- joined %>%
      filter(
        database %in% input$db,
        celltype %in% input$celltypes,
        str_detect(
          str_to_lower(as.character(feature_label)),
          "interferon|ifng|ifna|ifnb|type[-_ ]+(i|ii|1|2)[-_ ]+interferon"
        )
      ) %>%
      mutate(
        selected_program = biological_program %in% input$programs,
        selected_direction = case_when(
          input$direction == "Both" ~ TRUE,
          input$direction == "Up after treatment" ~ treatment_NES > 0,
          input$direction == "Down after treatment" ~ treatment_NES < 0,
          TRUE ~ TRUE
        ),
        selected_for_analysis = selected_program & selected_direction
      )

    recovery_keys <- recovery_direction_set() %>%
      distinct(pathway_row_key) %>%
      mutate(in_recovery = TRUE)

    step1_keys <- report_step1() %>%
      distinct(pathway_row_key) %>%
      mutate(in_step1 = TRUE)

    step2_keys <- report_step2() %>%
      distinct(pathway_row_key) %>%
      mutate(in_step2 = TRUE)

    step3_keys <- report_step3() %>%
      distinct(pathway_row_key) %>%
      mutate(in_step3 = TRUE)

    ifn_base %>%
      left_join(recovery_keys, by = "pathway_row_key") %>%
      left_join(step1_keys, by = "pathway_row_key") %>%
      left_join(step2_keys, by = "pathway_row_key") %>%
      left_join(step3_keys, by = "pathway_row_key") %>%
      mutate(
        across(
          c(in_recovery, in_step1, in_step2, in_step3),
          ~coalesce(.x, FALSE)
        ),
        classification_stage = case_when(
          !selected_for_analysis ~ "Outside current program/direction selection",
          in_recovery ~ "Removed as recovery",
          in_step3 ~ "Step 3 retained",
          in_step2 ~ "Step 2 retained; fails Step 3",
          in_step1 ~ "Step 1 retained; fails Step 2",
          is.na(treatment_FDR) ~ "Treatment comparison missing",
          treatment_FDR >= input$treat_fdr ~ "Not treatment-significant",
          TRUE ~ "Treatment-significant; not retained downstream"
        ),
        ifn_specificity = case_when(
          biological_program == "IFN-gamma / type II interferon signaling" ~
            "Type II / IFN-gamma",
          biological_program == "IFN-alpha/beta / type I interferon signaling" ~
            "Type I / IFN-alpha-beta",
          biological_program == "Broad / shared interferon signaling" ~
            "Broad/shared interferon parent",
          biological_program == "Antiviral response" ~
            "Antiviral / interferon-related",
          biological_program == "Other cytokine / JAK-STAT signaling" ~
            "Other cytokine/JAK-STAT",
          TRUE ~ "Other classification"
        )
      ) %>%
      select(
        celltype,
        database,
        feature_id,
        clean_pathway,
        biological_program,
        ifn_specificity,
        classification_stage,
        treatment_NES,
        treatment_FDR,
        disease_NES,
        disease_FDR,
        in_recovery,
        in_step1,
        in_step2,
        in_step3
      ) %>%
      distinct() %>%
      arrange(
        factor(
          biological_program,
          levels = c(
            "Broad / shared interferon signaling",
            "IFN-alpha/beta / type I interferon signaling",
            "IFN-gamma / type II interferon signaling",
            "Antiviral response",
            "Other cytokine / JAK-STAT signaling"
          )
        ),
        celltype,
        database,
        clean_pathway
      )
  })

  output$ifn_classification_qc_table <- renderDT({
    x <- interferon_classification_qc()

    validate(
      need(
        nrow(x) > 0,
        "No interferon-labeled pathways were found for the selected cell types/databases."
      )
    )

    standard_dt(x, page_length = 50)
  })

  report_step3 <- reactive({
    x <- report_step2() %>%
      filter(
        !is.na(abs_treatment_NES),
        abs_treatment_NES >= input$treatment_nes_min
      )

    if (isTRUE(input$use_disease_nes)) {
      x <- x %>%
        filter(
          !is.na(abs_disease_NES),
          abs_disease_NES <= input$disease_nes_max
        )
    }

    x
  })

  output$download_html_report <- downloadHandler(
    filename = function() {
      paste0(
        "PBMC_pathway_report_",
        format(Sys.Date(), "%Y%m%d"),
        ".html"
      )
    },

    content = function(file) {
      withProgress(message = "Generating PBMC HTML report...", value = 0, {

        req(length(input$db) > 0)
        req(length(input$celltypes) > 0)
        req(length(input$programs) > 0)

        incProgress(0.10, detail = "Collecting current analysis results...")

        report_data <- list(
          generated_at = Sys.time(),
          analysis_flow_gv = analysis_flow_gv,
          settings = list(
            treat_fdr = input$treat_fdr,
            disease_fdr_min = input$disease_fdr_min,
            exclude_recovery = input$exclude_recovery,
            treatment_nes_min = input$treatment_nes_min,
            use_disease_nes = input$use_disease_nes,
            disease_nes_max = input$disease_nes_max,
            db = input$db,
            celltypes = input$celltypes,
            programs = input$programs,
            direction = input$direction
          ),
          recovery = recovery_direction_set(),
          interferon_classification_qc = interferon_classification_qc(),
          step_contribution_blocks = step_contribution_blocks(),
          step1 = report_step1(),
          step2 = report_step2(),
          step3 = report_step3(),
          funnel = analysis_funnel()
        )

        incProgress(0.25, detail = "Creating temporary report files...")

        timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
        temp_report_dir <- file.path(
          tempdir(),
          paste0("pbmc_pathway_report_", timestamp)
        )
        dir.create(temp_report_dir, recursive = TRUE, showWarnings = FALSE)

        temp_qmd <- file.path(temp_report_dir, "PBMC_pathway_report.qmd")
        temp_rds <- file.path(temp_report_dir, "PBMC_pathway_report_data.rds")
        temp_html <- file.path(temp_report_dir, "PBMC_pathway_report.html")
        render_stdout <- file.path(temp_report_dir, "render_stdout.txt")
        render_stderr <- file.path(temp_report_dir, "render_stderr.txt")

        # The QMD is generated directly from this app; no external template file
        # is required in getwd().
        writeLines(generate_pbmc_qmd(), temp_qmd, useBytes = TRUE)
        saveRDS(report_data, temp_rds)

        if (!file.exists(temp_qmd) || !file.exists(temp_rds)) {
          stop("Failed to create temporary report inputs.")
        }

        incProgress(0.45, detail = "Checking Quarto...")

        if (!requireNamespace("quarto", quietly = TRUE)) {
          stop(
            "The R package 'quarto' is required. Install it with ",
            "install.packages('quarto')."
          )
        }

        if (!quarto::quarto_available()) {
          stop(
            "Quarto itself is not available on this computer. ",
            "Install Quarto, then restart R/RStudio."
          )
        }

        incProgress(0.60, detail = "Rendering self-contained HTML...")

        # Render in a clean temporary directory so the report is independent
        # of the user's working directory / OneDrive path.
        old_wd <- getwd()
        on.exit(setwd(old_wd), add = TRUE)
        setwd(temp_report_dir)

        render_error <- NULL

        tryCatch({
          if (requireNamespace("callr", quietly = TRUE)) {
            callr::r(
              function(qmd_file, rds_file, html_file) {
                if (!requireNamespace("quarto", quietly = TRUE)) {
                  stop("R package 'quarto' is not installed in the render process.")
                }

                quarto::quarto_render(
                  input = qmd_file,
                  output_file = basename(html_file),
                  execute_params = list(
                    report_rds = rds_file
                  ),
                  quiet = FALSE
                )
              },
              args = list(
                qmd_file = temp_qmd,
                rds_file = temp_rds,
                html_file = temp_html
              ),
              stdout = render_stdout,
              stderr = render_stderr
            )
          } else {
            quarto::quarto_render(
              input = temp_qmd,
              output_file = basename(temp_html),
              execute_params = list(
                report_rds = temp_rds
              ),
              quiet = FALSE
            )
          }
        }, error = function(e) {
          render_error <<- e
        })

        if (!is.null(render_error)) {
          stderr_text <- if (file.exists(render_stderr)) {
            paste(readLines(render_stderr, warn = FALSE), collapse = "\n")
          } else {
            ""
          }

          stdout_text <- if (file.exists(render_stdout)) {
            paste(readLines(render_stdout, warn = FALSE), collapse = "\n")
          } else {
            ""
          }

          stop(
            "HTML report rendering failed: ", render_error$message,
            if (nzchar(stderr_text)) paste0("\n\nRender stderr:\n", stderr_text) else "",
            if (nzchar(stdout_text)) paste0("\n\nRender stdout:\n", stdout_text) else ""
          )
        }

        incProgress(0.90, detail = "Preparing download...")

        if (!file.exists(temp_html)) {
          stop(
            "Quarto finished but the expected HTML file was not created: ",
            temp_html
          )
        }

        ok <- file.copy(temp_html, file, overwrite = TRUE)
        if (!isTRUE(ok)) {
          stop("Could not copy the rendered HTML to the download location.")
        }

        incProgress(1.00, detail = "Done")
      })
    },
    contentType = "text/html"
  )


  output$filter_text <- renderText({
    stage_note <- switch(
      analysis_stage(),
      step1 = paste0(
        "STEP 1 — RECOVERY PRE-FILTER\n",
        "Treatment FDR < ", input$treat_fdr, "\n",
        "Remove recovery when Disease NES and Treatment NES have opposite signs.\n",
        "Disease FDR is NOT used in Step 1.\n",
        "NES magnitude thresholds are NOT used in Step 1.\n"
      ),
      step2 = paste0(
        "STEP 2 — FDR CANDIDATE SET\n",
        "Start from Step 1 survivors.\n",
        "Require Treatment FDR < ", input$treat_fdr,
        " and Disease FDR >= ", input$disease_fdr_min, ".\n",
        "NES magnitude thresholds are NOT used in Step 2.\n"
      ),
      step3 = paste0(
        "STEP 3 — NES PRIORITIZATION\n",
        "Start from Step 2 candidates.\n",
        "Require |Treatment NES| >= ", input$treatment_nes_min, ".\n",
        if (isTRUE(input$use_disease_nes))
          paste0("Require |Disease NES| <= ", input$disease_nes_max, ".\n")
        else
          "No maximum |Disease NES| threshold is currently applied.\n"
      )
    )

    paste0(
      stage_note,
      "\nCurrent settings\n",
      "----------------\n",
      "Treatment FDR cutoff: ", input$treat_fdr, "\n",
      "Disease FDR minimum: ", input$disease_fdr_min, "\n",
      "Recovery removal enabled: ", input$exclude_recovery, "\n",
      "Recovery rule: opposite Disease/Treatment NES signs\n",
      "Step 3 minimum |Treatment NES|: ", input$treatment_nes_min, "\n",
      "Step 3 Disease NES guardrail enabled: ", input$use_disease_nes, "\n",
      if (isTRUE(input$use_disease_nes))
        paste0("Step 3 maximum |Disease NES|: ", input$disease_nes_max, "\n")
      else "",
      "Databases: ", paste(input$db, collapse = ", "), "\n",
      "Treatment direction: ", input$direction, "\n",
      "\nRecovery-direction rows identified: ", nrow(recovery_direction_set()), "\n",
      "Treatment-significant rows after recovery removal: ", nrow(recovery_cleaned_treatment_universe()), "\n",
      "Rows in current displayed result: ", nrow(filtered()), "\n",
      "Cell types remaining: ", n_distinct(filtered()$celltype), "\n",
      "Biological programs remaining: ", n_distinct(filtered()$biological_program), "\n"
    )
  })

}

shinyApp(ui, server)
