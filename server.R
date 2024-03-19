library(shiny)
library(shinyjs)
library(shinycssloaders)
library(DT)
library(tidyverse)
library(Matrix)
library(patchwork)
library(CCPlotR)
library(rvest)
library(BiocManager)
library(ggplot2)
library(viridis)
options(repos = BiocManager::repositories())
options(shiny.maxRequestSize = 100*1024^2)

# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {
  
  help <- readRDS("www/helptable.rds") # Read in dataframe of issues for help section
  
  genes_ids <- readRDS('www/genes_ids')
  
  scrape_gene_info <- function(gene_name) {
    gene_id = genes_ids[gene_name]
    url <- paste0('https://www.ncbi.nlm.nih.gov/gene/', gene_id)
    page <- rvest::read_html(url) %>% html_elements("dt:contains('Summary') + dd")
    summary <- html_text(page)
    return(c(summary,url))
  }
  
  data <- reactive({
    req(input$upload)
    inFile <- input$upload
    if (is.null(inFile)) { return(NULL) }
    fileExt <- tools::file_ext(inFile$name)
    datafile <- NULL
    if (fileExt == "rds") {
      datafile <- readRDS(inFile$datapath)
    } else {
      stop("Unsupported file type")
    }
  })
  sortedbycell <- reactive ({
    req(data(), input$cellchoice)
    data()[which(data()$source == input$cellchoice | data()$target == input$cellchoice), ]
  })
  
  output$interaction <- renderUI ({
    req(input$cellchoice, sortedbycell())
    selectInput('interaction', label = 'Select Interaction', choices = sort(unique(paste0(sortedbycell()$ligand, '|', sortedbycell()$receptor))))
  })
  
  output$celltype <- renderUI ({
    req(input$cellchoice, sortedbycell())
    selectInput('celltype', label = 'Select Cell Type', choices = sort(unique(c(sortedbycell()$source, '|', sortedbycell()$target))))
  })
  
  output$gene <- renderUI ({
    req(input$cellchoice, sortedbycell())
    selectInput('gene', label = 'Select Gene', choices = sort(unique(c(sortedbycell()$ligand, '|', sortedbycell()$receptor))), multiple = T)
  })
  
  
  exp_mtr <- reactive({
    req(input$expmtr)
    inFile <- input$expmtr
    if (is.null(inFile)) {
      return(NULL)
    }
    datafile <- readRDS(inFile$datapath)
  })
  
  meta <- reactive({
    req(input$meta)
    inFile <- input$meta
    if (is.null(inFile)) {
      return(NULL)
    }
    datafile <- readRDS(inFile$datapath)
  })
  
  cols <- reactive ({
    req(input$cols)
    inFile <- input$cols
    if(is.null(inFile)) {
      return(NULL)
    }
    datafile <- readRDS(inFile$datapath)
  })
  
  observe({
    updateSelectInput(session, "cellchoice", choices = sort(unique(c(data()$source, data()$target))))
  })
  
  output$dynamicCheckboxGroup <- renderUI({
    req(sortedbycell())
    df <- sortedbycell()  #
    colNames <- names(df) 
    checkboxGroupInput("show_cols", "Columns to display",
                       choices = colNames,
                       selected = colNames[1:4])
  })
  
  selected_columns <- reactive({
    req(input$show_cols) 
    input$show_cols
  })
  
  output$table <- DT::renderDT({
    req(input$upload)
    cols_to_show <- selected_columns()
    df <- sortedbycell()
    datatable(
      df,
      extensions = 'Buttons',
      options = list(
        autoWidth = TRUE, 
        scrollX = T, 
        buttons = c('csv', 'excel'), 
        dom = 'Bfrtip', 
        pageLength = 20,
        columnDefs = list(
          list(targets = which(names(df) %in% cols_to_show), visible = TRUE), # Show only selected columns
          list(targets = which(!names(df) %in% cols_to_show), visible = FALSE)
        ),
        initComplete = JS(
          "function(settings, json) {",
          "  var table = this.api();",
          "  table.on('click', 'td', function() {",
          "    var colIdx = table.cell(this).index().column;",
          "    if (table.column(colIdx).header().textContent === 'ligand' || table.column(colIdx).header().textContent === 'receptor') {",
          "      var geneName = table.cell(this).data();",
          "      Shiny.setInputValue('selected_gene', geneName);",
          "    } else if (table.column(colIdx).header().textContent === 'source' || table.column(colIdx).header().textContent === 'target') {",
          "      var cellName = table.cell(this).data();",
          "      Shiny.setInputValue('selected_cell', cellName);",
          "      }",
          "  });",
          "}")),
      filter = list(position = 'top', clear = TRUE),
      class = "display")
  })
  observeEvent(input$selected_gene, {
    gene_name <- input$selected_gene
    scraped_info <- scrape_gene_info(gene_name)
    gene_info <- scraped_info[1]
    title =  tags$a(href = scraped_info[2], target = "_blank", paste0("From NCBI: ", gene_name))
    showModal(modalDialog(
      title = title,
      gene_info,
      easyClose = TRUE,
      footer = NULL
    ))
  })

 
  output$int_plot <- renderPlot({
    lig <- str_extract(input$interaction, '[^|]+')
    rec <- str_extract(input$interaction, '[^|]+$')
    plot_df <- ungroup(sortedbycell()) %>% pivot_longer(cols = starts_with('agg'), names_to = 'tp', values_to = 'score', names_prefix = 'aggregate_rank_', values_drop_na = T) %>%
      mutate(tp = factor(tp, levels = c('Healthy', 'Diagnosis', 'Relapse', 'Post-treatment'), labels = c('Healthy', 'Diagnosis', 'Relapse', 'Post-Treatment')),
             from_cell = ifelse(source == input$cellchoice, paste('From', input$cellchoice), paste('To', input$cellchoice))) %>%
      filter((source == input$cellchoice | target == input$cellchoice) & ligand == lig & receptor == rec)
    if(input$plot_type_int == 'Heatmap'){
      print(cc_heatmap(plot_df %>% filter(ligand == lig, receptor == rec), option = 'B') + 
              facet_grid(from_cell ~ tp, scales = 'free_x', switch = 'y', space = 'free_x') +
              theme(strip.placement = 'outside', legend.key.height = unit(4.5, 'lines')))
    } 
    if(input$plot_type_int == 'Connections'){
      cell_cols <- cols()$celltype
      h_plot <- cc_sigmoid(plot_df %>% filter(ligand == lig, receptor == rec, tp == 'Healthy'), colours = cell_cols) + 
        labs(title = 'Healthy') +
        theme(plot.title = element_text(hjust = 0.5, face = 'bold'))
      d_plot <- cc_sigmoid(plot_df %>% filter(ligand == lig, receptor == rec, tp == 'Diagnosis'), colours = cell_cols) + 
        labs(title = 'Diagnosis') +
        theme(plot.title = element_text(hjust = 0.5, face = 'bold'))
      r_plot <- cc_sigmoid(plot_df %>% filter(ligand == lig, receptor == rec, tp == 'Relapse'), colours = cell_cols) + 
        labs(title = 'Relapse') +
        theme(plot.title = element_text(hjust = 0.5, face = 'bold'))
      pt_plot <- cc_sigmoid(plot_df %>% filter(ligand == lig, receptor == rec, tp == 'Post-Treatment'), colours = cell_cols) + 
        labs(title = 'Post-Treatment') +
        theme(plot.title = element_text(hjust = 0.5, face = 'bold'))
      print(h_plot + d_plot + r_plot + pt_plot)
    }
    if(input$plot_type_int == 'Chord diagram'){
      cell_cols <- cols()$celltype
      par(oma = c(4,1,1,1), mfrow = c(2, 2), mar = c(2, 2, 1, 1))
        if(nrow(plot_df %>% filter(ligand == lig, receptor == rec, tp == 'Healthy')) > 0){
          cc_circos(ungroup(plot_df) %>% filter(ligand == lig, receptor == rec, tp == 'Healthy'), cell_cols = cell_cols, option = 'B', cex = 1, show_legend = F, scale = T)
          title('Healthy')}
        if(nrow(plot_df %>% filter(ligand == lig, receptor == rec, tp == 'Diagnosis')) > 0){
          cc_circos(ungroup(plot_df) %>% filter(ligand == lig, receptor == rec, tp == 'Diagnosis'), cell_cols = cell_cols, option = 'B', cex = 1, show_legend = F, scale = T) 
          title('Diagnosis')}
        if(nrow(plot_df %>% filter(ligand == lig, receptor == rec, tp == 'Relapse')) > 0){
          cc_circos(ungroup(plot_df) %>% filter(ligand == lig, receptor == rec, tp == 'Relapse'), cell_cols = cell_cols, option = 'B', cex = 1, show_legend = F, scale = T) 
          title('Relapse')}
        if(nrow(plot_df %>% filter(ligand == lig, receptor == rec, tp == 'Post-Treatment')) > 0){
          cc_circos(ungroup(plot_df) %>% filter(ligand == lig, receptor == rec, tp == 'Post-Treatment'), cell_cols = cell_cols, option = 'B', cex = 1, show_legend = F, scale = T) 
          title('Post-Treatment')}
      par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
      plot(1, type = "n", axes=FALSE, xlab="", ylab="")
      legend(x = "bottom", horiz = F,
             legend = unique(c((plot_df %>% filter(ligand == lig, receptor == rec) %>% pull(source)), (plot_df %>% filter(ligand == lig, receptor == rec) %>% pull(target)))),
             title = "Cell type",
             pch = 15,
             ncol = ceiling(length(unique(c((plot_df %>% filter(ligand == lig, receptor == rec) %>% pull(source)), (plot_df %>% filter(ligand == lig, receptor == rec) %>% pull(target)))))/2),
             text.width = max(sapply(unique(c((plot_df %>% filter(ligand == lig, receptor == rec) %>% pull(source)), (plot_df %>% filter(ligand == lig, receptor == rec) %>% pull(target)))), strwidth)),
             xpd = TRUE,
             col = cell_cols[unique(c((plot_df %>% filter(ligand == lig, receptor == rec) %>% pull(source)), (plot_df %>% filter(ligand == lig, receptor == rec) %>% pull(target))))])
    }
    if (input$plot_type_int == "Violin plot") {
      req(exp_mtr())
      req(meta())
      expdf <- exp_mtr() 
      meta <- meta()
      cell_cols <- cols()$celltype
      exp_df <- cbind(meta, data.frame(lig = expdf[,lig], rec = expdf[,rec]))
      p1 <- ggplot(exp_df, aes(x = cell_type, y = lig, fill = timepoint)) +
        geom_violin(show.legend = F, scale = 'width', col = 'black', draw_quantiles = 0.5) +
        scale_fill_manual(values = cols()$timepoint) +
        scale_x_discrete(limits = names(cell_cols)) +
        labs(y = lig) +
        theme_classic(base_size = 18) +
        theme(axis.text = element_text(colour = 'black'),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.line.x = element_blank(),
              axis.title.x = element_blank(),
              legend.position = 'none')
      p2 <- ggplot(exp_df, aes(x = cell_type, y = rec)) +
        geom_violin(aes(fill = timepoint), scale = 'width', draw_quantiles = 0.5) +
        scale_fill_manual(values = cols()$timepoint, name= 'Timepoint') +
        scale_x_discrete(limits = names(cell_cols)) +
        guides(colour = guide_legend(override.aes = list(size = 3))) +
        labs(y = rec) +
        theme_classic(base_size = 18) +
        theme(axis.text = element_text(colour = 'black'),
              axis.text.x  = element_text(angle = 90, hjust = 1, vjust = 0.5),
              axis.title.x = element_blank(),
              legend.position = 'bottom')
      print(p1/p2)
    }
  })
  
  output$cell_plot <- renderPlot({
    plot_df <- sortedbycell() %>% pivot_longer(cols = starts_with('agg'), names_to = 'tp', values_to = 'score', names_prefix = 'aggregate_rank_', values_drop_na = T) %>%
      mutate(tp = factor(tp, levels = c('Healthy', 'Diagnosis', 'Relapse', 'Post-treatment'), labels = c('Healthy', 'Diagnosis', 'Relapse', 'Post-Treatment')))
    
    if(input$plot_type_cell == 'Heatmap'){
      p1 <- cc_heatmap(ungroup(plot_df) %>% filter(((source == input$cellchoice & target == input$celltype) | (source == input$celltype & target == input$cellchoice)), tp == 'Healthy'), option = 'B', n_top_ints = input$n_ints_cell) + 
        scale_fill_viridis_c(option = 'C', na.value = 'black', direction = 1, limits=plot_df %>% 
                               filter((source == input$cellchoice & target == input$celltype) | (source == input$celltype & target == input$cellchoice)) %>%
                               group_by(tp) %>% slice_max(order_by = score, n = input$n_ints_cell) %>% 
                               pull(score) %>% range()) +
        labs(title = 'Healthy') +
        theme(plot.title = element_text(hjust = 0.5),
              legend.key.height = unit(0.3, 'inches'))
      
      p2 <- cc_heatmap(ungroup(plot_df) %>% filter(((source == input$cellchoice & target == input$celltype) | (source == input$celltype & target == input$cellchoice)), tp == 'Diagnosis'), option = 'B', n_top_ints = input$n_ints_cell) + 
        scale_fill_viridis_c(option = 'C', na.value = 'black', direction = 1, limits=plot_df %>% 
                               filter((source == input$cellchoice & target == input$celltype) | (source == input$celltype & target == input$cellchoice)) %>%
                               group_by(tp) %>% slice_max(order_by = score, n = input$n_ints_cell) %>% 
                               pull(score) %>% range()) +
        labs(title = 'Diagnosis') +
        theme(plot.title = element_text(hjust = 0.5),
              legend.key.height = unit(0.3, 'inches'))
      p3 <- cc_heatmap(ungroup(plot_df) %>% filter(((source == input$cellchoice & target == input$celltype) | (source == input$celltype & target == input$cellchoice)), tp == 'Relapse'), option = 'B', n_top_ints = input$n_ints_cell) + 
        scale_fill_viridis_c(option = 'C', na.value = 'black', direction = 1, limits=plot_df %>% 
                               filter((source == input$cellchoice & target == input$celltype) | (source == input$celltype & target == input$cellchoice)) %>%                               group_by(tp) %>% slice_max(order_by = score, n = input$n_ints_cell) %>%
                               pull(score) %>% range()) +
        labs(title = 'Relapse') +
        theme(plot.title = element_text(hjust = 0.5),
              legend.key.height = unit(0.3, 'inches'))
      p4 <- cc_heatmap(ungroup(plot_df) %>% filter(((source == input$cellchoice & target == input$celltype) | (source == input$celltype & target == input$cellchoice)), tp == 'Post-Treatment'), option = 'B', n_top_ints = input$n_ints_cell) + 
        scale_fill_viridis_c(option = 'C', na.value = 'black', direction = 1, limits=plot_df %>% 
                               filter((source == input$cellchoice & target == input$celltype) | (source == input$celltype & target == input$cellchoice)) %>%
                               group_by(tp) %>% slice_max(order_by = score, n = input$n_ints_cell) %>% 
                               pull(score) %>% range()) +
        labs(title = 'Post-Treatment') +
        theme(plot.title = element_text(hjust = 0.5),
              legend.key.height = unit(0.3, 'inches'))     
      
      print(p1+p2+p3+p4 + plot_layout(guides = 'collect'))
    }
    if(input$plot_type_cell == 'Connections'){
      cell_cols <- cols()$celltype
      h_plot <- cc_sigmoid(plot_df %>% filter(((source == input$cellchoice & target == input$celltype) | (source == input$celltype & target == input$cellchoice)), tp == 'Healthy'), colours = cell_cols, n_top_ints = input$n_ints_cell) +
        labs(title = 'Healthy') +
        theme(plot.title = element_text(hjust = 0.5, face = 'bold'))
      d_plot <- cc_sigmoid(plot_df %>% filter(((source == input$cellchoice & target == input$celltype) | (source == input$celltype & target == input$cellchoice)), tp == 'Diagnosis'), colours = cell_cols, n_top_ints = input$n_ints_cell) +
        labs(title = 'Diagnosis') +
        theme(plot.title = element_text(hjust = 0.5, face = 'bold'))
      r_plot <- cc_sigmoid(plot_df %>% filter(((source == input$cellchoice & target == input$celltype) | (source == input$celltype & target == input$cellchoice)), tp == 'Relapse'), colours = cell_cols, n_top_ints = input$n_ints_cell) +
        labs(title = 'Relapse') +
        theme(plot.title = element_text(hjust = 0.5, face = 'bold'))
      pt_plot <- cc_sigmoid(plot_df %>% filter(((source == input$cellchoice & target == input$celltype) | (source == input$celltype & target == input$cellchoice)), tp == 'Post-Treatment'), colours = cell_cols, n_top_ints = input$n_ints_cell) +
        labs(title = 'Post-Treatment') +
        theme(plot.title = element_text(hjust = 0.5, face = 'bold'))
      print(h_plot + d_plot)
    }
    if(input$plot_type_cell == 'Chord diagram'){
      cell_cols <- cols()$celltype
      par(oma = c(4,1,1,1), mfrow = c(2, 2), mar = c(2, 2, 1, 1))
      healthy_df <- plot_df %>% 
        filter(((source == input$cellchoice & target == input$celltype) | 
                  (source == input$celltype & target == input$cellchoice)) & 
                 tp == 'Healthy')
      diagnosis_df <- plot_df %>% 
        filter(((source == input$cellchoice & target == input$celltype) | 
                  (source == input$celltype & target == input$cellchoice)) & 
                 tp == 'Diagnosis')
      relapse_df <- plot_df %>% 
        filter(((source == input$cellchoice & target == input$celltype) | 
                  (source == input$celltype & target == input$cellchoice)) & 
                 tp == 'Relapse')
      post_t_df <- plot_df %>% 
        filter(((source == input$cellchoice & target == input$celltype) | 
                  (source == input$celltype & target == input$cellchoice)) & 
                 tp == 'Post-Treatment')
      if(nrow(healthy_df) > 0){
        cc_circos(ungroup(healthy_df), cell_cols = cell_cols, option = 'B', cex = 0.8, show_legend = F, scale = T, n_top_ints = input$n_ints_cell)
        title('Healthy')}
      if(nrow(diagnosis_df) > 0){
        cc_circos(ungroup(diagnosis_df), cell_cols = cell_cols, option = 'B', cex = 0.8, show_legend = F, scale = T, n_top_ints = input$n_ints_cell)
        title('Diagnosis')}
      if(nrow(relapse_df) > 0){
        cc_circos(ungroup(relapse_df), cell_cols = cell_cols, option = 'B', cex = 0.8, show_legend = F, scale = T, n_top_ints = input$n_ints_cell)
        title('Relapse')}
      if(nrow(post_t_df) > 0){
        cc_circos(ungroup(post_t_df), cell_cols = cell_cols, option = 'B', cex = 0.8, show_legend = F, scale = T, n_top_ints = input$n_ints_cell)
        title('Post-Treatment')}
      par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
      plot(1, type = "n", axes=FALSE, xlab="", ylab="")
      unique_cells <- unique(c(healthy_df$source, healthy_df$target, diagnosis_df$source, diagnosis_df$target, relapse_df$source, relapse_df$target, post_t_df$source, post_t_df$target))
      legend(x = "bottom", legend = unique_cells, title = "Cell type", pch = 15, ncol = ceiling(length(unique_cells)/2), text.width = max(sapply(unique_cells, strwidth)), xpd = TRUE, col = cell_cols[unique_cells])
    }
  })
  
  output$h_netplot <- renderPlot({
    plot_df <- sortedbycell() %>% pivot_longer(cols = starts_with('agg'), names_to = 'tp', values_to = 'score', names_prefix = 'aggregate_rank_', values_drop_na = T) %>%
      mutate(tp = factor(tp, levels = c('Healthy', 'Diagnosis', 'Relapse', 'Post-treatment'), labels = c('Healthy', 'Diagnosis', 'Relapse', 'Post-Treatment')))
    cell_cols <- cols()$celltype
    if(input$plot_type_cell == 'Network diagram'){
      cc_network(plot_df %>% filter((source == input$cellchoice & target == input$celltype), tp == 'Healthy'), colours = cell_cols, n_top_ints = input$n_ints_cell, option = 'B', node_size = 2.2) +
        labs(title = 'Healthy') +
        theme(plot.title = element_text(hjust = 0.5, face = 'bold'))
    }
  })
  
  output$d_netplot <- renderPlot({
    plot_df <- sortedbycell() %>% pivot_longer(cols = starts_with('agg'), names_to = 'tp', values_to = 'score', names_prefix = 'aggregate_rank_', values_drop_na = T) %>%
      mutate(tp = factor(tp, levels = c('Healthy', 'Diagnosis', 'Relapse', 'Post-treatment'), labels = c('Healthy', 'Diagnosis', 'Relapse', 'Post-Treatment')))
    cell_cols <- cols()$celltype
    if(input$plot_type_cell == 'Network diagram'){
      cc_network(plot_df %>% filter((source == input$cellchoice & target == input$celltype), tp == 'Diagnosis'), colours = cell_cols, n_top_ints = input$n_ints_cell, option = 'B', node_size = 2.2) +
        labs(title = 'Diagnosis') +
        theme(plot.title = element_text(hjust = 0.5, face = 'bold'))
    }
  })
  
  output$r_netplot <- renderPlot({
    plot_df <- sortedbycell() %>% pivot_longer(cols = starts_with('agg'), names_to = 'tp', values_to = 'score', names_prefix = 'aggregate_rank_', values_drop_na = T) %>%
      mutate(tp = factor(tp, levels = c('Healthy', 'Diagnosis', 'Relapse', 'Post-treatment'), labels = c('Healthy', 'Diagnosis', 'Relapse', 'Post-Treatment')))
    cell_cols <- cols()$celltype
    if(input$plot_type_cell == 'Network diagram'){
      cc_network(plot_df %>% filter((source == input$cellchoice & target == input$celltype), tp == 'Relapse'), colours = cell_cols, n_top_ints = input$n_ints_cell, option = 'B', node_size = 2.2) +
        labs(title = 'Relapse') +
        theme(plot.title = element_text(hjust = 0.5, face = 'bold'))
    }
  })
  
  output$pt_netplot <- renderPlot({
    plot_df <- sortedbycell() %>% pivot_longer(cols = starts_with('agg'), names_to = 'tp', values_to = 'score', names_prefix = 'aggregate_rank_', values_drop_na = T) %>%
      mutate(tp = factor(tp, levels = c('Healthy', 'Diagnosis', 'Relapse', 'Post-treatment'), labels = c('Healthy', 'Diagnosis', 'Relapse', 'Post-Treatment')))
    cell_cols <- cols()$celltype
    if(input$plot_type_cell == 'Network diagram'){
      cc_network(plot_df %>% filter((source == input$cellchoice & target == input$celltype), tp == 'Post-Treatment'), colours = cell_cols, n_top_ints = input$n_ints_cell, option = 'B', node_size = 2.2) +
        labs(title = 'Post-Treatment') +
        theme(plot.title = element_text(hjust = 0.5, face = 'bold'))
    }
  })
  
  output$networkplot <- renderUI({
    if(input$plot_type_cell == 'Network diagram'){
      tagList(
      fluidRow(
        column(6, plotOutput("h_netplot", height = "400px")),
        column(6, plotOutput("d_netplot", height = "400px"))
      ),
      fluidRow(
        column(6, plotOutput("r_netplot", height = "400px")),
        column(6, plotOutput("pt_netplot", height = "400px"))
      )
      )
    }
  })
  
  
  output$gene_plot <- renderPlot({
    genes <- input$gene
    plot_df <- sortedbycell() %>% pivot_longer(cols = starts_with('agg'), names_to = 'tp', values_to = 'score', names_prefix = 'aggregate_rank_', values_drop_na = T) %>%
      mutate(tp = factor(tp, levels = c('Healthy', 'Diagnosis', 'Relapse', 'Post-treatment'), labels = c('Healthy', 'Diagnosis', 'Relapse', 'Post-Treatment'))) %>%
      filter(((source == input$cellchoice | target == input$cellchoice)) & 
               (ligand %in% genes | receptor %in% genes))
    if(input$plot_type_gene == 'Heatmap'){
      p1 <- cc_heatmap(ungroup(plot_df) %>% filter((ligand %in% genes | receptor %in% genes), tp == 'Healthy'), option = 'B', n_top_ints = input$n_ints_gene) + 
        scale_fill_viridis_c(option = 'C', na.value = 'black', direction = 1, limits=plot_df %>% 
                               filter(ligand %in% genes | receptor %in% genes) %>% 
                               group_by(tp) %>% slice_max(order_by = score, n = input$n_ints_gene) %>% 
                               pull(score) %>% range()) +
        labs(title = 'Healthy') +
        theme(plot.title = element_text(hjust = 0.5),
              legend.key.height = unit(0.3, 'inches'))
      
      p2 <- cc_heatmap(ungroup(plot_df) %>% filter((ligand %in% genes | receptor %in% genes), tp == 'Diagnosis'), option = 'B', n_top_ints = input$n_ints_gene) + 
        scale_fill_viridis_c(option = 'C', na.value = 'black', direction = 1, limits=plot_df %>% 
                               filter(ligand %in% genes | receptor %in% genes) %>% 
                               group_by(tp) %>% slice_max(order_by = score, n = input$n_ints_gene) %>% 
                               pull(score) %>% range()) +
        labs(title = 'Diagnosis') +
        theme(plot.title = element_text(hjust = 0.5),
              legend.key.height = unit(0.3, 'inches'))
      p3 <- cc_heatmap(ungroup(plot_df) %>% filter((ligand %in% genes | receptor %in% genes), tp == 'Relapse'), option = 'B', n_top_ints = input$n_ints_gene) + 
        scale_fill_viridis_c(option = 'C', na.value = 'black', direction = 1, limits=plot_df %>% 
                               filter(ligand %in% genes | receptor %in% genes) %>% 
                               group_by(tp) %>% slice_max(order_by = score, n = input$n_ints_gene) %>% 
                               pull(score) %>% range()) +
        labs(title = 'Relapse') +
        theme(plot.title = element_text(hjust = 0.5),
              legend.key.height = unit(0.3, 'inches'))
      p4 <- cc_heatmap(ungroup(plot_df) %>% filter((ligand %in% genes | receptor %in% genes), tp == 'Post-Treatment'), option = 'B', n_top_ints = input$n_ints_gene) + 
        scale_fill_viridis_c(option = 'C', na.value = 'black', direction = 1, limits=plot_df %>% 
                               filter(ligand %in% genes | receptor %in% genes) %>% 
                               group_by(tp) %>% slice_max(order_by = score, n = input$n_ints_gene) %>% 
                               pull(score) %>% range()) +
        labs(title = 'Post-treatment') +
        theme(plot.title = element_text(hjust = 0.5),
              legend.key.height = unit(0.3, 'inches'))
      print(p1+p2+p3+p4 + plot_layout(guides = 'collect'))
    }
    if(input$plot_type_gene == 'Connections'){
      cell_cols <- cols()$celltype
      h_plot <- cc_sigmoid(plot_df %>% filter((ligand %in% genes | receptor %in% genes) & tp == 'Healthy'), colours = cell_cols, n_top_ints = input$n_ints_gene) +
        labs(title = 'Healthy') +
        theme(plot.title = element_text(hjust = 0.5, face = 'bold'))
      d_plot <- cc_sigmoid(plot_df %>% filter((ligand %in% genes | receptor %in% genes) & tp == 'Diagnosis'), colours = cell_cols, n_top_ints = input$n_ints_gene) +
        labs(title = 'Diagnosis') +
        theme(plot.title = element_text(hjust = 0.5, face = 'bold'))
      r_plot <- cc_sigmoid(plot_df %>% filter((ligand %in% genes | receptor %in% genes) & tp == 'Relapse'), colours = cell_cols, n_top_ints = input$n_ints_gene) +
        labs(title = 'Relapse') +
        theme(plot.title = element_text(hjust = 0.5, face = 'bold'))
      pt_plot <- cc_sigmoid(plot_df %>% filter((ligand %in% genes | receptor %in% genes) & tp == 'Post-Treatment'), colours = cell_cols, n_top_ints = input$n_ints_gene) +
        labs(title = 'Post-Treatment') +
        theme(plot.title = element_text(hjust = 0.5, face = 'bold'))
      print(h_plot + d_plot + r_plot + pt_plot)
    }
    if(input$plot_type_gene == 'Chord diagram'){
      cell_cols <- cols()$celltype
      par(oma = c(4,1,1,1), mfrow = c(2, 2), mar = c(2, 2, 1, 1))
      if(nrow(plot_df %>% filter((ligand %in% genes | receptor %in% genes) & tp == 'Healthy')) > 0){
        cc_circos(ungroup(plot_df) %>% filter((ligand %in% genes | receptor %in% genes) & tp == 'Healthy'), cell_cols = cell_cols, option = 'B', cex = 0.8, show_legend = F, scale = T, n_top_ints = input$n_ints_gene)
        title('Healthy')}
      if(nrow(plot_df %>% filter((ligand %in% genes | receptor %in% genes) & tp == 'Diagnosis')) > 0){
        cc_circos(ungroup(plot_df) %>% filter((ligand %in% genes | receptor %in% genes) & tp == 'Diagnosis'), cell_cols = cell_cols, option = 'B', cex = 0.8, show_legend = F, scale = T, n_top_ints = input$n_ints_gene)
        title('Diagnosis')}
      if(nrow(plot_df %>% filter((ligand %in% genes | receptor %in% genes) & tp == 'Relapse')) > 0){
        cc_circos(ungroup(plot_df) %>% filter((ligand %in% genes | receptor %in% genes) & tp == 'Relapse'), cell_cols = cell_cols, option = 'B', cex = 0.8, show_legend = F, scale = T, n_top_ints = input$n_ints_gene)
        title('Relapse')}
      if(nrow(plot_df %>% filter((ligand %in% genes | receptor %in% genes) & tp == 'Post-Treatment')) > 0){
        cc_circos(ungroup(plot_df) %>% filter((ligand %in% genes | receptor %in% genes) & tp == 'Post-Treatment'), cell_cols = cell_cols, option = 'B', cex = 0.8, show_legend = F, scale = T, n_top_ints = input$n_ints_gene)
        title('Post-Treatment')}
      par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
      plot(1, type = "n", axes=FALSE, xlab="", ylab="")
      legend(x = "bottom", horiz = F,
             legend = unique(c((plot_df %>% filter(ligand %in% genes | receptor %in% genes) %>% pull(source)), (plot_df %>% filter(ligand %in% genes | receptor %in% genes) %>% pull(target)))),
             title = "Cell type",
             pch = 15,
             ncol = ceiling(length(unique(c((plot_df %>% filter(ligand %in% genes | receptor %in% genes) %>% pull(source)), (plot_df %>% filter(ligand %in% genes | receptor %in% genes) %>% pull(target)))))/2),
             text.width = max(sapply(unique(c((plot_df %>% filter(ligand %in% genes | receptor %in% genes) %>% pull(source)), (plot_df %>% filter(ligand %in% genes | receptor %in% genes) %>% pull(target)))), strwidth)),
             xpd = TRUE,
             col = cell_cols[unique(c((plot_df %>% filter(ligand %in% genes | receptor %in% genes) %>% pull(source)), (plot_df %>% filter(ligand %in% genes | receptor %in% genes) %>% pull(target))))])
    }
    if(input$plot_type_gene == 'Violin plot'){
      req(exp_mtr())
      req(meta())
      expdf <- exp_mtr()
      meta <- meta()
      cell_cols <- cols()$celltype
      exp_df <- cbind(meta, (as.data.frame(as.matrix(expdf[, input$gene])) %>% setNames(input$gene))) %>% pivot_longer(input$gene, names_to = 'gene', values_to = 'value')
      ggplot(exp_df, aes(x = cell_type, y = value, fill = timepoint)) +
        geom_violin(show.legend = T, scale = 'width', col = 'black', draw_quantiles = 0.5) +
        scale_fill_manual(values = cols()$timepoint, name = 'Timepoint') +
        scale_x_discrete(limits = names(cell_cols)) +
        labs(y = 'Normalised expression', x = NULL) +
        facet_grid(gene~., switch = 'y') +
        theme_classic(base_size = 18) +
        theme(axis.text = element_text(colour = 'black'),
              axis.text.x = element_text(hjust=1, angle = 90, vjust = 0.5),
              strip.placement = 'outside',
              legend.position = 'bottom')
    }
  })
  
  combined_data <- reactive ({
    req(data())
    data() %>%
      select(source, target, ligand, receptor, timepoints) %>%
      pivot_longer(cols = c(source, target), names_to = "column", values_to = "cell_type") })
  
  
  
  data_with_indicators <- reactive ({
    req(combined_data())
    combined_data() %>%
      mutate(
        HealthyPresent = as.integer(str_detect(timepoints, "Healthy")),
        DiagnosisPresent = as.integer(str_detect(timepoints, "Diagnosis")),
        RelapsePresent = as.integer(str_detect(timepoints, "Relapse")),
        PostTreatmentPresent = as.integer(str_detect(timepoints, "Post-treatment"))
      )
  })
  
  counts_by_condition <- reactive({
    req(data_with_indicators())
    total_obs <- nrow(data_with_indicators())
    
    bind_rows(
      data_with_indicators() %>%
        filter(HealthyPresent == 1) %>%
        group_by(cell_type) %>%
        summarise(Count = n(), .groups = "drop") %>%
        group_by(Condition = "Healthy") %>%
        mutate(Proportion = Count / sum(Count)),
      data_with_indicators() %>%
        filter(DiagnosisPresent == 1) %>%
        group_by(cell_type) %>%
        summarise(Count = n(), .groups = "drop") %>%
        group_by(Condition = "Diagnosis") %>%
        mutate(Proportion = Count / sum(Count)),
      data_with_indicators() %>%
        filter(RelapsePresent == 1) %>%
        group_by(cell_type) %>%
        summarise(Count = n(), .groups = "drop") %>%
        group_by(Condition = "Relapse") %>%
        mutate(Proportion = Count / sum(Count)),
      data_with_indicators() %>%
        filter(PostTreatmentPresent == 1) %>%
        group_by(cell_type) %>%
        summarise(Count = n(), .groups = "drop") %>%
        group_by(Condition = "Post-Treatment") %>%
        mutate(Proportion = Count / sum(Count))
    )
  })
  
  output$summary_plot_prop <- renderPlotly({
    if(input$summary_plot == 'Cell Proportion') {
      plot_ly(counts_by_condition(), x = ~cell_type, y = ~Proportion, color = ~Condition, type = 'bar') %>%
        layout(
          title = "Cell Proportion Plot by Condition",
          xaxis = list(title = "Cell Type"),
          yaxis = list(title = "Proportion"),
          hovermode = 'closest'
        )
    }
  })

  
output$summary_plot_int <- renderPlot ({
  if(input$summary_plot == 'Top Interactions by Cell') {
    req(sortedbycell())


    plot_df <- sortedbycell() %>%
      pivot_longer(cols = starts_with('agg'), names_to = 'tp', values_to = 'score', names_prefix = 'aggregate_rank_', values_drop_na = TRUE) %>%
      mutate(tp = factor(tp, levels = c('Healthy', 'Diagnosis', 'Relapse', 'Post-treatment'), labels = c('Healthy', 'Diagnosis', 'Relapse', 'Post-Treatment')),
             Pair = paste(ligand, receptor, sep = " - "))

    unique_pairs <- unique(plot_df$Pair)

    num_colors <- length(unique_pairs)
    random_colors <- sample(viridis(num_colors), num_colors)

    colour_palette <- setNames(random_colors, unique_pairs)

    plots_list <- list()


    for (tp in levels(plot_df$tp)) {

      top_pairs <- plot_df %>%
        filter(tp == !!tp) %>%
        count(Pair, sort = TRUE) %>%
        head(10)  # Get the top 10 pairs


      p <- ggplot(top_pairs, aes(x = reorder(Pair, n), y = n, fill = Pair)) +
        geom_bar(stat = "identity") +
        coord_flip() +
        scale_fill_manual(values = colour_palette) +
        labs(title = paste("Top 10 Ligand-Receptor Pairs in", input$cellchoice, "Cells - ", tp), x = "Ligand - Receptor Pair", y = "Frequency")

      plots_list[[tp]] <- p
    }

    combined_plot <- wrap_plots(plots_list, ncol = 2)
    print(combined_plot)
  }

})

observeEvent(input$summary_plot, {
  if (input$summary_plot == "Change in Interactions") {
    output$change_int <- renderUI({
      selectInput('interactions', 'Select Interaction', 
                  choices = sort(unique(paste0(sortedbycell()$ligand, '|', sortedbycell()$receptor))))
    })
  } else {
    output$dynamic_ui <- renderUI({}) # Clear the UI if other options are selected
  }
})

output$summary_plot_change <- renderPlotly ({
  req(input$interactions)
  plot_df <- sortedbycell() %>% pivot_longer(cols = starts_with('agg'), names_to = 'tp', values_to = 'score', names_prefix = 'aggregate_rank_', values_drop_na = T) %>%
    mutate(tp = factor(tp, levels = c('Healthy', 'Diagnosis', 'Relapse', 'Post-treatment'), labels = c('Healthy', 'Diagnosis', 'Relapse', 'Post-Treatment')),
           Pair = paste(ligand, receptor, sep = "|"))
  print(head(plot_df$Pair))
  
  selectedPair <- input$interactions
  print(selectedPair)
  
  filtered <- plot_df %>%
    filter(Pair == selectedPair)
 print(filtered$Pair)
  
  totalpertimepoint <- plot_df %>%
    group_by(tp) %>%
    summarise(total=n())
  
  countpertimepoint <- filtered %>%
    group_by(tp) %>%
    summarise(count = n())
  
  proportionData <- inner_join(countpertimepoint, totalpertimepoint, by = "tp") %>%
    mutate(proportion = count / total)
print(head(proportionData))
  
  gg <- ggplot(proportionData, aes(x = tp, y = proportion, group = 1)) +
    geom_line() +
    geom_point() +
    theme_classic() +
    theme(plot.title = element_text(face = "bold", size = 14),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 10)) +
    labs(title = paste("Proportion of", selectedPair, "Interaction per Timepoint"),
         x = "Timepoint",
         y = "Proportion")
  
  # Convert to plotly object
  p <- ggplotly(gg)
  
  # Display the plotly plot
  
})

  
  selected_title <- reactive({
    help_topic <- input$help
    index <- which(help$title == help_topic)
    help$issue[index]
  })
  
  selected_help <- reactive({
    help_topic <- input$help
    index <- which(help$title == help_topic)
    help$comment[index]
  })
  
  output$help_issue <- renderText({
    selected_title()
  })
  
  output$help_comment <- renderUI({
    eval(parse(text =selected_help()))
  })
  
  
  session$onSessionEnded(stopApp)
})