library(shiny)
library(shinyjs)
library(DT)
library(shinycssloaders)
library(BiocManager)
library(ggplot2)
library(plotly)
options(repos = BiocManager::repositories())

help <- readRDS("www/helptable.rds") # Read in list of issues for help section
issues <- help$title

shinyUI(fluidPage(
  
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "custom.css")
  ),
  useShinyjs(),
  titlePanel(tags$div(class='flex-container',
                      tags$div(class="flex-child", img(src ='interactome.png', width = '20%')),
                      tags$div(class="flex-child", tags$h1("Cell-Cell Interactome"),
                               tags$p(tags$b('Explore cell-cell interactions between different cell types from scRNA-seq data', style = "font-size:18px;")),
                               tags$p(tags$b('Citation:'), 'Cell-cell interactome of the hematopoietic niche and its changes in acute myeloid leukemia.', tags$i("Ennis S et. al.,", 'iScience, 2023. DOI:', tags$a('10.1016/j.isci.2023.106943', href = 'https://doi.org/10.1016/j.isci.2023.106943')), style = "font-size:16px;"),
                               tags$p(tags$b('Created by: '), tags$a(' Sarah Ennis', href='mailto:s.ennis6@universityofgalway.ie', .noWS='outside'), ',', tags$a(' Micheál Ó Dálaigh', href='mailto:m.odalaigh1@universityofgalway.ie', .noWS='outside'), ',', tags$a(' Jacopo Umberto Verga', href='mailto:j.verga1@universityofgalway.ie', .noWS='outside'), ',',tags$a(' Pilib Ó Broin', href='mailto:pilib.obroin@universityofgalway.ie', .noWS='outside'),',', tags$a(' Eva Szegezdi', href='mailto:eva.szegezdi@universityofgalway.ie', .noWS='outside'), ' and', tags$a(' Siobhan Griffin', href='mailto:s.griffin25@universityofgalway.ie', .noWS='outside'), style = "font-size:16px;", .noWS = c("after-begin", "before-end")))),
             windowTitle = 'HSC Interactome'),
  
  fluidRow( 
    column(4, fileInput("upload", "Interactome Data", accept = c("rds"))),
    column(4, fileInput("cols", "Colour Palette", accept = c("rds")))),
  selectInput("cellchoice", label = 'Select cell type', choices = sort(unique(c(data()$source, data()$target)))),
  
  tabsetPanel(id='panels',
              tabPanel('Table',
                       div(class = 'mainDiv',
                           sidebarLayout(
                             sidebarPanel(
                               width = 3,
                               uiOutput("dynamicCheckboxGroup")),
                               #checkboxGroupInput("show_cols", "Columns to display", choices = NULL)), #(choiceValues = 1:ncol(data())),
                             #  choiceNames = colnames(data())),
                             mainPanel(
                               withSpinner(DT::DTOutput('table'), color="#9FDC93", size = 0.5))
                           ))),
              tabPanel('Plots',
                       div(class = 'mainDiv',
                           sidebarLayout(
                             sidebarPanel(
                               width = 3,
                               fileInput("expmtr", "Expression Matrix (For Violin Plots only)", accept = c("rds")), 
                               fileInput("meta", "Meta Data (For Violin Plots only)", accept = c("rds")),
                               tabsetPanel(id="tabs",
                                           tabPanel('Interactions',
                                                    br(),
                                                    uiOutput('interaction'),
                                                    hr(),
                                                    radioButtons('plot_type_int', 'Select plot type', choices = c('Heatmap', 'Connections', 'Chord diagram', 'Violin plot'), selected = 'Heatmap'),
                                                    hr()
                                           ),
                                           tabPanel('Cell types',
                                                    br(),
                                                    uiOutput('celltype'),
                                                    hr(),
                                                    radioButtons('plot_type_cell', 'Select plot type', choices = c('Heatmap', 'Connections', 'Chord diagram', 'Network diagram'), selected = 'Heatmap'),
                                            
                                                    hr(),
                                                    sliderInput('n_ints_cell', 'Max no. of interactions to plot', min = 5, max = 40, value = 20, step = 1)
                                           ),
                                           tabPanel('Genes',
                                                    br(),
                                                    uiOutput('gene'), 
                                                    hr(),
                                                    radioButtons('plot_type_gene', 'Select plot type', choices = c('Heatmap', 'Connections', 'Chord diagram', 'Violin plot'), selected = 'Heatmap'),
                                                    hr(),
                                                    sliderInput('n_ints_gene', 'Max no. of interactions to plot', min = 5, max = 40, value = 20, step = 1)
                                           ))
                             ),
                             mainPanel(
                               conditionalPanel(condition="input.tabs == 'Interactions'",
                                                withSpinner(plotOutput('int_plot', width = '100%', height = "600px"), color="#9FDC93", size = 0.5)
                               ),
                               conditionalPanel(condition="input.tabs == 'Cell types' && (input.plot_type_cell == 'Heatmap' || input.plot_type_cell == 'Connections' || input.plot_type_cell == 'Chord diagram')",
                                                withSpinner(plotOutput('cell_plot', width = '100%', height = "600px"), color="#9FDC93", size = 0.5)
                               ),                
                               conditionalPanel(condition = "input.tabs == 'Cell types' && input.plot_type_cell == 'Network diagram'",
                                                withSpinner(uiOutput("networkplot"), color="#9FDC93", size = 0.5)
                                                ),
                               conditionalPanel(condition="input.tabs == 'Genes'",
                                                withSpinner(plotOutput('gene_plot', width = '100%', height = "600px"), color="#9FDC93", size = 0.5)
                               )
                             )))
              ),
              tabPanel('Summary',
                       div(class = 'mainDiv',
                           sidebarLayout(
                             sidebarPanel(
                               width = 3,
                               tabsetPanel(id="tabs2",
                                           tabPanel("Summary Type",
                                                    br(),
                                                    radioButtons('summary_plot', 'Select summary type', choices = c('Cell Proportion', 'Top Interactions by Cell', 'Change in Interactions'), selected = 'Cell Proportion'),
                                                    uiOutput("change_int")
                                           )
                               )),
                             mainPanel(
                               conditionalPanel(
                                 condition = "input.summary_plot == 'Cell Proportion'",
                                 withSpinner(plotlyOutput('summary_plot_prop', width = '100%', height = "600px"), color="#9FDC93", size = 0.5)
                               ),
                               conditionalPanel(
                                 condition = "input.summary_plot == 'Top Interactions by Cell'",
                                 withSpinner(plotOutput('summary_plot_int', width = '100%', height = "600px"), color="#9FDC93", size = 0.5)
                               ),
                               conditionalPanel(
                                 condition = "input.summary_plot == 'Change in Interactions'",
                                 withSpinner(plotlyOutput('summary_plot_change', width = '100%', height = "600px"), color="#9FDC93", size = 0.5)
                               )        
                               
                             )))
              ),
              
              tabPanel('Help',
                       div(class = 'helpDiv',
                           tags$li('Welcome to the help section of the Cell-Cell Interactome R Shiny app'),
                           tags$li('Please see the drop-down list below for common questions regarding the dataset and how to navigate the app'),
                           tags$li('If your question is not answered here, please open an issue on our GitHub issues page ', tags$a('here', href = 'https://github.com/SzegezdiLab/HSC_Interactome/issues')),
                           
                           
                           selectInput("help", label = h3("Please select an issue:", style = "font-weight: bold;"),
                                       choices = issues),
                           
                           # Render help section based on issue chosen
                           tags$h3(textOutput("help_issue"), style = "font-weight: bold;"),
                           tags$div(uiOutput("help_comment"))
                       )
              ))))