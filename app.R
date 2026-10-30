# setup and loading packages ----------------------------------------------

library("shiny")
library("ggplot2")
library("DESeq2")
library("mosdef")
library("dplyr")
library("GeneTonic")
library("shinythemes")
library("shinycssloaders")
library("plotly")
library("patchwork")

# dataset config 

default_genes <- c("FDXR", "MDM2", "CDKN1A", "PVT1", "PMAIP1")

# dataset 1: (X-ray only)
dds_phybion <- readRDS("data/dds_xray_old.RDS")
dds_phybion <- estimateSizeFactors(dds_phybion)

# dataset 2: neutron vs. X-ray
dds_nx <- readRDS("data/dds_neutron_xray.RDS") 
colnames(colData(dds_nx))[6] <- "Donors"
dds_nx$Dose <- factor(dds_nx$Dose, levels = c("0Gy", "0.1Gy", "0.25Gy", "0.5Gy", "1Gy", "2Gy", "4Gy"))
design(dds_nx) <- ~Donors + group
dds_nx <- estimateSizeFactors(dds_nx)
dds_nx$Time_point <- factor(dds_nx$Time_point, levels = c("6h", "24h"))

nx_time_choices <- levels(droplevels(dds_nx$Time_point))

nx_radiation_col <- "Radiation_type"  
nx_time_col      <- "Time_point"      
nx_dose_col      <- "Dose"           

dds_list_phybion <- list(
  `2hrs` = dds_phybion[, dds_phybion$Time_point == "2hrs"],
  `6hrs` = dds_phybion[, dds_phybion$Time_point == "6hrs"]
)

# annotation per dataset (in case gene sets / IDs differ)
anno_list <- list(
  phybion = mosdef::get_annotation_orgdb(dds_phybion, "org.Hs.eg.db", id_type = "ENSEMBL"),
  nx      = mosdef::get_annotation_orgdb(dds_nx, "org.Hs.eg.db", id_type = "ENSEMBL")
)

# choices for the neutron/X-ray dropdowns (taken from the data)
nx_radiation_choices <- sort(unique(as.character(colData(dds_nx)[[nx_radiation_col]])))
# global color palettes: one fixed color per level, across both datasets ------
color_vars <- c("Sex", "Donors", "Dose")

all_levels <- function(var) {
  u <- unique(unlist(lapply(list(dds_phybion, dds_nx), function(d) {
    if (var %in% colnames(colData(d))) as.character(colData(d)[[var]])
  })))
  if (var == "Dose") {
    num <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", u)))
    u <- u[order(!is.na(num), num)]
  } else {
    u <- sort(u)
  }
  u
}

make_palette <- function(lv) {
  base <- RColorBrewer::brewer.pal(8, "Set2")
  cols <- if (length(lv) <= 8) base[seq_along(lv)] else colorRampPalette(base)(length(lv))
  setNames(cols, lv)
}

color_palettes <- lapply(setNames(color_vars, color_vars), function(v) make_palette(all_levels(v)))

# fixed sex colors: male = blue, female = pink
sex_lv <- all_levels("Sex")
color_palettes$Sex <- setNames(
  ifelse(grepl("^f", sex_lv, ignore.case = TRUE), "#E75480", "#4A90D9"),
  sex_lv
)

# defining helper functions -----------------------------------------------

plot_gene_modified <- function(dds_list, genes, intgroup = "condition", assay = "counts",
                               annotation_obj = NULL, normalized = TRUE, transform = TRUE,
                               color_by = "Donors", plot_titles = names(dds_list)) {
  
  if (!is.null(plot_titles) && length(plot_titles) != length(dds_list)) {
    stop("The length of `plot_titles` must match the number of elements in `dds_list`.")
  }
  
  df_list <- lapply(seq_along(dds_list), function(i) {
    dds <- dds_list[[i]]
    
    data_list <- lapply(genes, function(gene) {
      df <- mosdef::get_expr_values(de_container = dds, gene = gene, intgroup = intgroup,
                                    assay = assay, normalized = normalized)
      df$sample_id <- rownames(df)
      df$gene_id <- gene
      df$source <- plot_titles[i]
      
      if (color_by %in% colnames(colData(dds))) {
        df[[color_by]] <- colData(dds)[, color_by]
      } else {
        df[[color_by]] <- NA
      }
      
      return(df)
    })
    
    bind_rows(data_list)
  })
  
  df <- bind_rows(df_list)
  df$source <- factor(df$source, levels = unique(plot_titles))
  df$plotby <- interaction(df[[intgroup]])
  
  if (!is.null(annotation_obj)) {
    df$gene_name <- annotation_obj$gene_name[match(df$gene_id, annotation_obj$gene_id)]
  } else {
    df$gene_name <- df$gene_id
  }
  
  plots <- lapply(unique(df$gene_name), function(gene_name) {
    gene_df <- df[df$gene_name == gene_name, ]
    
    ggplot(gene_df, aes(x = .data$plotby, y = .data$exp_value, color = .data[[color_by]],
                        group = .data$source,
                        text = paste0("Sample: ", sample_id,
                                      "<br>Value: ", round(exp_value, 2),
                                      "<br>Time: ", source,
                                      "<br>Dose: ", .data$plotby))) +
      geom_point(position = position_jitter(width = 0.2), shape = 16, size = 3, alpha = 0.8) +
      stat_summary(fun = mean, geom = "line", aes(group = source), color = "grey70", linewidth = 1) +
      facet_wrap(~source) +
      theme_minimal(base_size = 14) +
      labs(title = gene_name, x = "", y = ifelse(transform, "Normalized counts (log10)", "Normalized counts"),
           color = color_by) +
      scale_color_manual(values = color_palettes[[color_by]], na.value = "grey50") +
      {
        if (transform) scale_y_log10() else scale_y_continuous()
      }
  })
  
  return(plots)
}

# subset the neutron/X-ray dds by radiation type and timepoint(s);
# one list element per timepoint so the existing faceting is reused
make_nx_dds_list <- function(dds, radiation, timepoint) {
  cd <- colData(dds)
  tps <- if (identical(timepoint, "All")) nx_time_choices else timepoint
  
  out <- lapply(tps, function(tp) {
    keep <- as.character(cd[[nx_radiation_col]]) == radiation &
      as.character(cd[[nx_time_col]]) == tp
    dds[, keep]
  })
  names(out) <- tps
  Filter(function(x) ncol(x) > 0, out)
}


# here comes the app definition ------------------------------------------------

## ui definition ---------------------------------------------------------------

phybion_ui <- fluidPage(
  theme = shinytheme("flatly"),
  title = "phybion_explorer: Exploring Gene Dynamics in Response to Radiation",
  titlePanel(div(
    h2("📈 Exploring Gene Dynamics in Response to Radiation"),
    style = "margin-bottom: 20px"
  )),
  
  sidebarLayout(
    sidebarPanel(
      wellPanel(
        selectInput("dataset", label = "Dataset",
                    choices = c("D1_Xray(2h&6h)" = "phybion",
                                "D2_Neutron_X-ray(6h&24h)"            = "nx")),
        
        # only shown for the neutron/X-ray dataset
        conditionalPanel(
          condition = "input.dataset == 'nx'",
          selectInput("radiation_type", label = "Radiation type",
                      choices = nx_radiation_choices),
          selectInput("timepoint", label = "Timepoint",
                      choices = c("All", nx_time_choices))
        ),
        
        selectizeInput("gene_name", label = "Select gene(s)", choices = NULL, multiple = TRUE),
        selectInput("color_by", label = "Color by", choices = c("Sex", "Donors", "Dose")),
        downloadButton("download_plot", "Download Plot"),
        uiOutput("genecard_links")
      )
    ),
    
    mainPanel(
      h4("Here you can visualize normalized expression across doses and timepoints."),
      withSpinner(uiOutput("dynamic_plot_ui"), type = 4)
    )
  )
)

## server definition -----------------------------------------------------------
phybion_server <- function(input, output, session) {
  
  # annotation for the currently selected dataset
  current_anno <- reactive({
    anno_list[[input$dataset]]
  })
  
  # repopulate gene choices whenever the dataset changes
  observeEvent(input$dataset, {
    anno <- current_anno()
    updateSelectizeInput(
      session = session,
      inputId = "gene_name",
      choices = c(Choose = "", anno$gene_name),
      selected = intersect(default_genes, anno$gene_name),
      server = TRUE
    )
  })
  
  # dds list for the selected dataset / radiation type / timepoint
  current_dds_list <- reactive({
    if (input$dataset == "phybion") {
      dds_list_phybion
    } else {
      req(input$radiation_type, input$timepoint)
      lst <- make_nx_dds_list(dds_nx, input$radiation_type, input$timepoint)
      validate(need(length(lst) > 0,
                    "No samples for this radiation type / timepoint combination."))
      lst
    }
  })
  
  # Store ggplot objects
  gene_plot_data <- reactive({
    req(input$gene_name)
    anno <- current_anno()
    gene_ids <- anno$gene_id[match(input$gene_name, anno$gene_name)]
    
    plot_gene_modified(
      dds_list = current_dds_list(),
      genes = gene_ids,
      intgroup = "Dose",
      annotation_obj = anno,
      color_by = input$color_by
    )
  })
  
  gene_plots <- reactive({
    lapply(gene_plot_data(), ggplotly, tooltip = "text")
  })
  
  output$dynamic_plot_ui <- renderUI({
    req(input$gene_name)
    plot_ids <- seq_along(input$gene_name)
    
    plot_outputs <- lapply(plot_ids, function(i) {
      plotlyOutput(outputId = paste0("gene_plot_", i), height = "500px")
    })
    
    do.call(tagList, plot_outputs)
  })
  
  observe({
    req(input$gene_name)
    plots <- gene_plots()
    
    lapply(seq_along(plots), function(i) {
      output[[paste0("gene_plot_", i)]] <- renderPlotly({
        plots[[i]]
      })
    })
  })
  
  output$download_plot <- downloadHandler(
    filename = function() {
      prefix <- if (input$dataset == "nx") paste0(input$radiation_type, "_", input$timepoint, "_") else ""
      paste0(prefix, paste(input$gene_name, collapse = "_"), "_expression_plot.png")
    },
    content = function(file) {
      ggsave(file, plot = wrap_plots(gene_plot_data(), ncol = 1),
             device = "png", width = 12, height = 4 * length(input$gene_name))
    }
  )
  
  output$genecard_links <- renderUI({
    req(input$gene_name)
    links <- lapply(input$gene_name, function(gene) {
      gene_link <- paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=", gene)
      tags$a(href = gene_link, target = "_blank", paste("GeneCards:", gene))
    })
    tagList(tags$h5("GeneCards Links:"), tags$ul(lapply(links, tags$li)))
  })
}

# run the application ----------------------------------------------------------
shinyApp(ui = phybion_ui, server = phybion_server)
