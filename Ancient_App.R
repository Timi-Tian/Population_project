rm(list=ls())
library(shiny)
library(dplyr)
library(DT)
library(readr)

country_list = c("Default","Argentina", "Armenia", "Austria", "Azerbaijan", "Bahamas", 
                 "Belgium", "Bolivia", "Botswana", "Cameroon", "Canary Islands",
                 "Channel Islands", "China", "Croatia", "Cuba", "Czech Republic",
                 "DR Congo", "Denmark",  "Estonia", "Faroes", "France",       
                 "Germany", "Greenland", "Guadeloupe", "Hungary", "Iceland",    
                 "Iran", "Ireland", "Isle of Man", "Israel", "Italy",       
                 "Japan", "Jordan", "Kazakhstan", "Kenya", "Kyrgyzstan",   
                 "Latvia", "Malawi", "Mongolia", "Nepal", "Netherlands",   
                 "Norway", "Peru", "Poland",  "Portugal", "Puerto Rico", 
                 "Romania", "Russia", "South Africa", "Spain", "St. Lucia",    
                 "Sweden", "Switzerland", "Syria", "Taiwan", "Tanzania",    
                 "Turkey", "Turkmenistan", "USA", "Ukraine", "United Kingdom")

hap_ID <- read_csv("/Users/tianmi/Desktop/unique_haplogroups.csv")

ui <- fluidPage(
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "https://cdnjs.cloudflare.com/ajax/libs/font-awesome/5.15.1/css/all.min.css")
  ),
  tags$head(
    tags$style(HTML("
      .dt-head-center { text-align: center; }
      th { background-color: lightblue; color: black; }
    "))
  ),
  # App title
  titlePanel("Haplogroup Path"),
  # Sidebar layout with a input
  sidebarLayout(
    position = "left",
    # Sidebar panel for inputs
    sidebarPanel(
      # Input: Selector for choosing haplogroups
      selectInput(inputId = "Haplogroups",
                  choices = hap_ID$Unique_Haplogroups,
                  label = "Select a haplogroup ID:"),
      selectInput(inputId = "Country",
                  choices = country_list,
                  label = "Select a Country:",
                  selected = "Default" ) # default selection
    ),
    # Main panel containing elements that get created in the server function
    mainPanel(
      htmlOutput(outputId = "final_path"),
      dataTableOutput(outputId = "final_table")
    )
  )
)


# Define server logic
server <- function(input, output) {
  result <- reactive({
    # code to generate the final table here
    source('~/Desktop/Path.R', local = TRUE)
    data <- generate_final_table(input$Haplogroups, input$Country, merged_data)
    final_table <- data$table
    final_path <- data$path
    
    col_titles <- c("Haplogroup" = "Name of the ancestral haplogroup.",
                    "Age Estimate" = "The estimated time when the most recent ancestor of this lineage was born.",
                    "Archaeological Era" = "Global archaeological time period associated with the age estimate.",
                    "Time Passed" = "Time elapsed between this haplogroup and its ancestral haplogroup. A large number can suggest a small population size or a bottleneck, causing only one lineage to survive for a long time.",
                    "Immediate Descendants" = "Number of phylogenetic subclades. A large number indicates a rapid expansion event.",
                    "Tested Modern Desendants" = "The number of present-day DNA testers confirmed to belong to this haplogroup.")
    new_colnames <- sapply(names(col_titles), function(cn) {
      paste(cn, sprintf('<span title="%s">?</span>', col_titles[[cn]]))
    })
    names(final_table) <- new_colnames
    return(list(path = final_path, table = final_table))
  })

  output$final_path <- renderUI({
    final_path <- result()$path
    # Create a custom HTML output for the path
    HTML(paste("<p style='font-weight:bold; font-size:20px;'>Your detailed haplogroup path is:", final_path, "</p>"))
  })
  # Render the input haplogroup ID as output
  output$final_table <- renderDataTable({
    datatable(result()$table, escape = FALSE, options = list(
      autoWidth = TRUE,
      columnDefs = list(
        list(targets = "_all", className = 'dt-head-center')
      )
      )) %>%
      formatStyle(
        columns = names(result()$table),
        backgroundColor = styleEqual(levels = names(result()), values = rep('lightblue', length(names(result())))),
        color = 'black'
      )
  })
}
# Run the app
shinyApp(ui = ui, server = server)
