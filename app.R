library(shiny)
library(ggplot2)
library(dplyr)
library(stringr)
library(tibble)

ui <- fluidPage(
  titlePanel("Interactive QC Filtering for sc/snRNA-seq"),
  wellPanel(
    h4("📘 Usage Example"),
    p("This app expects a dataframe with the following structure:"),
    tableOutput("usage_example")
  ),
  sidebarLayout(
    sidebarPanel(
      radioButtons("data_source", "Choose Data Source:",
                   choices = c("Upload CSV", "Upload RDS")),
      conditionalPanel(
        condition = "input.data_source == 'Upload CSV'",
        fileInput("file_csv", "Upload CSV File")
      ),
      conditionalPanel(
        condition = "input.data_source == 'Upload RDS'",
        fileInput("file_rds", "Upload RDS File")
      ),
      hr(),
      selectInput("sample_mode", "Filter Mode:", choices = c("All Samples", "Single Sample")),
      uiOutput("sample_selector"),
      selectInput("source_filter", "Source Filter:", choices = c("raw", "filtered"), selected = "raw"),
      tags$div(
        style = "margin-top: 20px;",
        h4("Filtering Parameters"),
        wellPanel(
          h5("nCount_RNA"),
          fluidRow(
            column(6, numericInput("min_counts_num", "Min nCount_RNA:", value = 1000)),
            column(6, numericInput("max_counts_num", "Max nCount_RNA:", value = 10000))
          )
        ),
        wellPanel(
          h5("nFeature_RNA"),
          fluidRow(
            column(6, numericInput("min_features_num", "Min nFeature_RNA:", value = 500)),
            column(6, numericInput("max_features_num", "Max nFeature_RNA:", value = 5000))
          )
        ),
        wellPanel(
          h5("percent.mt"),
          fluidRow(
            column(6, numericInput("min_mito", "Min percent.mt:", value = 0)),
            column(6, numericInput("max_mito", "Max percent.mt:", value = 10))
          )
        )
      )
    ),
    mainPanel(
      h4(textOutput("cell_count")),
      fluidRow(
        column(4, h5("Remaining Cells per Sample"), tableOutput("sample_counts")),
        column(4, h5("Filtered-Out Cells per Sample"), tableOutput("filtered_counts")),
        column(4, h5("Percent Filtered-Out per Sample"), tableOutput("percent_filtered"))
      ),
      plotOutput("scatterPlot", height = "600px")
    )
  )
)

server <- function(input, output, session) {
  output$usage_example <- renderTable({
    tibble(
      CellID = c("TACGCTCCACCAATTG", "GCGGAAATCAGTCATG", "GGGATGATCCAGCAAT", "TATTCCATCGTACACA", "CTGTAGAGTTGTGGCC"),
      sample = factor(c("LGC-31758-sn", "LGC-31758-sn", "LGC-31663-sn", "HGH-31828-sn", "LGC-31758-sc")),
      source = c("raw", "raw", "filtered", "raw", "filtered"),
      nFeature_RNA = c(163, 135, 166, 905, 5568),
      nCount_RNA = c(465, 326, 222, 1918, 45036),
      percent.mt = c(0.430, 0.0000125, 1.35, 2.45, 0.966),
      sample_type = c("snRNA-seq", "snRNA-seq", "snRNA-seq", "snRNA-seq", "scRNA-seq")
    )
  })
  
  selected_df <- reactive({
    req(input$data_source)
    if (!is.null(input$file_csv)) {
      read.csv(input$file_csv$datapath)
    } else if (!is.null(input$file_rds)) {
      readRDS(input$file_rds$datapath)
    } else {
      NULL
    }
  })
  
  processed_df <- reactive({
    req(selected_df())
    df <- selected_df() %>%
      mutate(cell_id = row_number()) %>%
      mutate(sample_type = case_when(
        str_detect(sample, "-sc") ~ "scRNA-seq",
        str_detect(sample, "-sn") ~ "snRNA-seq",
        TRUE ~ "unknown"
      ))
    sample_order <- df %>%
      distinct(sample, sample_type) %>%
      arrange(factor(sample_type, levels = c("scRNA-seq", "snRNA-seq", "unknown"))) %>%
      pull(sample)
    df$sample <- factor(df$sample, levels = sample_order)
    df
  })
  
  output$sample_selector <- renderUI({
    req(processed_df())
    if (input$sample_mode == "Single Sample") {
      selectInput("selected_sample", "Choose Sample:", choices = unique(processed_df()$sample))
    }
  })
  
  filtered_data <- reactive({
    req(processed_df())
    df <- processed_df()
    if ("source" %in% colnames(df)) {
      df <- df %>% filter(source == input$source_filter)
    }
    df <- df %>%
      filter(
        nCount_RNA >= input$min_counts_num,
        nCount_RNA <= input$max_counts_num,
        nFeature_RNA >= input$min_features_num,
        nFeature_RNA <= input$max_features_num,
        percent.mt >= input$min_mito,
        percent.mt <= input$max_mito
      )
    if (input$sample_mode == "Single Sample") {
      req(input$selected_sample)
      df <- df %>% filter(sample == input$selected_sample)
    }
    df
  })
  
  output$sample_counts <- renderTable({
    req(filtered_data())
    df <- filtered_data() %>%
      group_by(sample) %>%
      summarise(Remaining_Cells = n(), .groups = "drop")
    total <- data.frame(sample = "Total", Remaining_Cells = sum(df$Remaining_Cells))
    bind_rows(total, df)
  })
  
  output$filtered_counts <- renderTable({
    req(processed_df(), filtered_data())
    df_all <- processed_df()
    if ("source" %in% colnames(df_all)) {
      df_all <- df_all %>% filter(source == input$source_filter)
    }
    df_filtered <- filtered_data()
    df_removed <- df_all %>% filter(!cell_id %in% df_filtered$cell_id)
    df <- df_removed %>%
      group_by(sample) %>%
      summarise(Filtered_Cells = n(), .groups = "drop")
    total <- data.frame(sample = "Total", Filtered_Cells = sum(df$Filtered_Cells))
    bind_rows(total, df)
  })
  
  output$percent_filtered <- renderTable({
    req(processed_df(), filtered_data())
    df_all <- processed_df()
    if ("source" %in% colnames(df_all)) {
      df_all <- df_all %>% filter(source == input$source_filter)
    }
    df_filtered <- filtered_data()
    df_summary <- df_all %>%
      group_by(sample) %>%
      summarise(Total_Cells = n(), .groups = "drop")
    df_filtered_summary <- df_filtered %>%
      group_by(sample) %>%
      summarise(Remaining_Cells = n(), .groups = "drop")
    df_percent <- left_join(df_summary, df_filtered_summary, by = "sample") %>%
      mutate(Remaining_Cells = ifelse(is.na(Remaining_Cells), 0, Remaining_Cells),
             Percent_Filtered = round(100 * (1 - Remaining_Cells / Total_Cells), 2)) %>%
      select(sample, Percent_Filtered)
    total_percent <- data.frame(
      sample = "Total",
      Percent_Filtered = round(mean(df_percent$Percent_Filtered, na.rm = TRUE), 2)
    )
    bind_rows(total_percent, df_percent)
  })
  
  output$scatterPlot <- renderPlot({
    req(processed_df())
    df_all <- processed_df()
    if ("source" %in% colnames(df_all)) {
      df_all <- df_all %>% filter(source == input$source_filter)
    }
    df_filtered <- filtered_data()
    df_all$filtered <- ifelse(df_all$cell_id %in% df_filtered$cell_id, "kept", "removed")
    if (input$sample_mode == "Single Sample") {
      req(input$selected_sample)
      df_all <- df_all %>% filter(sample == input$selected_sample)
    }
    p <- ggplot(df_all, aes(x = nCount_RNA, y = nFeature_RNA)) +
      geom_point(data = subset(df_all, filtered == "removed"), color = "grey80", alpha = 0.4, size = 0.6) +
      geom_point(data = subset(df_all, filtered == "kept"), aes(color = percent.mt), alpha = 0.6, size = 0.8) +
      scale_x_log10() +
      scale_y_log10() +
      scale_color_viridis_c(option = "plasma") +
      geom_vline(xintercept = c(input$min_counts_num, input$max_counts_num), linetype = "dashed", color = "red") +
      geom_hline(yintercept = c(input$min_features_num, input$max_features_num), linetype = "dashed", color = "blue") +
      theme_bw() +
      labs(
        x = "nCount_RNA (log scale)",
        y = "nFeature_RNA (log scale)",
        color = "Percent Mito"
      ) + geom_smooth()
    if (input$sample_mode == "All Samples") {
      p <- p + facet_wrap(~sample, scales = "free")
    }
    p
  })
}

shinyApp(ui = ui, server = server)
