library(ggpubr)
library(readxl)
library(caret)
library(dplyr)
library(tidyr)
library(plotly)
library(gt)
library(webshot2)
library(shinythemes)

# Define the UI
ui <- fluidPage(
  titlePanel("Bradford Analysis"),
  theme = shinytheme("cerulean"),
  sidebarLayout(
    sidebarPanel(
      fileInput("file", "Choose Excel File", accept = c(".xlsx")),
      sliderInput("micrograms", "Micrograms Sample", min = 0, max = 100, value = 60),
      textInput("dilution", "Dilution", value = "1.5"),
      actionButton("openLink", "Sample.xlsx"),
      actionButton("runButton", "Run Analysis"),
      downloadButton("downloadTable", "Download Table")
    ),
    mainPanel(
      plotOutput("brd_plot"),
      gt_output("table")
    )
  )
)

# Define the server
server <- function(input, output) {
  
  observeEvent(input$openLink, {
    browseURL("https://www.example.com")  # Replace with your desired web link
  })
  
  # Reactive value for micrograms_sample
  micrograms_sample <- reactive({
    input$micrograms
  })
  
  # Read the Excel file and generate brd_plot and table
  observeEvent(input$file, {
    bradford_raw <- read_excel(input$file$datapath, col_names = FALSE)
    bradford_train<-data.frame(absorbance=unlist(bradford_raw[1:9,2:3]), concentration=0)
    bradford_train[c(2,11), 2]<- 5
    bradford_train[c(3,12), 2]<- 10
    bradford_train[c(4,13), 2]<- 15
    bradford_train[c(5,14), 2]<- 20
    bradford_train[c(6,15), 2]<- 25
    bradford_train[c(7,16), 2]<- 30
    bradford_train[c(8,17), 2]<- 35
    bradford_train[c(9,18), 2]<- 40
    brd_line <- lm(concentration~absorbance, bradford_train)
    # ... Code for brd_plot ...
    # Load the necessary libraries
   # for stat_cor and stat_regline_equation
    
    # Set theme options for a professional look
    theme_set(theme_bw(base_size = 12, base_family = "Arial"))
    
    # Create the plot
    brd_plot <- ggplot(bradford_train, aes(x = concentration, y = absorbance)) +
      geom_smooth(method = "lm", formula = y ~ x, se = FALSE) +
      geom_point() +
      stat_cor(label.y = 0.9, label.x = 0.1) +
      stat_regline_equation(label.y = 0.95, label.x = 0.1) +
      labs(x = "Concentration", y = "Absorbance", title = "Bradford Training Data")
    
    # Customize the plot theme
    brd_plot +
      theme(
        plot.title = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        legend.position = "bottom"
      )
    
    
    
    # ... Code for final_table ...
    final_table <- bradford_raw[, 4:6]
    names(final_table) <- c("sample", "abs1", "abs2")
    final_table$absorbance <- rowMeans(final_table[, c("abs1", "abs2")], na.rm = TRUE)
    final_table <- final_table[, c("sample", "absorbance")]
    #`predicting concentration`
    final_table$concentrations<-predict.lm(brd_line, final_table)  # Replace with your code
    
    # Apply the dilutionobserveEvent(input$runButton, {
      micrograms_sample <- reactive({
        input$micrograms
      })
      
      dilution <- as.numeric(input$dilution)
  
    
    final_table$concentrations <- final_table$concentrations * dilution
    
    # Calculate volumes based on micrograms_sample
    final_table$volume_sample <- micrograms_sample() / final_table$concentrations
    final_table$laemli4x <- 0.25 * max(final_table$volume_sample) / 0.75
    final_table$final_volume <- 0.25 * max(final_table$volume_sample) / 0.75 + max(final_table$volume_sample)
    final_table$ripa <- final_table$final_volume - (final_table$volume_sample + final_table$laemli4x)
    
    # Rearrange the columns
    final_table <- final_table %>%
      mutate_if(is.numeric, round, digits = 3) %>%
      relocate(sample, absorbance, concentrations, laemli4x, ripa, volume_sample, final_volume)
  
    names(final_table)<-c("Sample","Absorbance","µg/µL protein","Laemli Solution","RIPA","Volume of Sample","Final Volume")
    total_laemli<-sum(final_table$`Laemli Solution`)
    # Create the gt table
    tbl <- gt(final_table)
   
    
    # ... Code for styling the table ...
    tbl <- tbl %>%
      tab_header(
        title = md("*Volumes Chart for Western Blots*"),
        subtitle = md("Everything is in **µL**")  # Change subtitle text color to blue
      ) %>%
      tab_style(
        style = cell_text(weight = "bold"),
        locations = cells_body(columns = 1)
      ) %>%
      tab_style(
        style = list(cell_fill(color = "#add8e6"),
                     cell_text(weight = "bold")),
        locations = cells_body(columns = 4:6))%>%
      tab_style(
        style = list(cell_fill(color = "#ffefb5"),
                     cell_text(weight ="bold")), 
        locations = cells_body(columns = 7))%>%
      tab_footnote(
        footnote=paste(round((total_laemli)/0.75), "µL of Laemli Solution will be needed. Add",round((total_laemli)/0.75)*0.1,"µL of Mercapto 10% into",round((total_laemli)/0.75)*0.9, "µL of Laemli 4X", sep=" "),
        locations = NULL,
        placement = "auto"
      )
    
    output$brd_plot <- renderPlot({
      brd_plot
    })
    
    observeEvent(input$runButton, {
      output$table <- render_gt({
       tbl
      })
    })
  })
  
  # Download the table as PNG
  # Download the table as PNG
  output$downloadTable <- downloadHandler(
    filename = function() {
      paste("table", Sys.Date(), ".png", sep = "_")
    },
    content = function(file) {
      # Save the table as a temporary HTML file
      tmp_file <- tempfile(fileext = ".html")
      gtsave(tbl, file = tmp_file)
      
      # Use webshot2 to capture a screenshot of the HTML file
      webshot2::webshot(url = tmp_file, file = file, delay = 2)
      
      # Delete the temporary HTML file
      file.remove(tmp_file)
    }
  )
  
}
# Run the application
shinyApp(ui = ui, server = server)
