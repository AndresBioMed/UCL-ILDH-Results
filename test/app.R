library(shiny)
library(ggplot2)
library(gt)
library(readxl)
library(webshot2)

# Define the UI
ui <- fluidPage(
  titlePanel("Bradford Analysis"),
  sidebarLayout(
    sidebarPanel(
      fileInput("file", "Choose Excel File", accept = c(".xlsx")),
      sliderInput("micrograms", "Micrograms Sample", min = 0, max = 100, value = 60),
      textInput("dilution", "Dilution", value = "1.5"),
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
    brd_plot <- ggplot(bradford_train, aes(x=concentration,y=absorbance))+geom_smooth(method = "lm", formula = y ~ x) +geom_point() + stat_cor(label.y = 0.9)+ 
      stat_regline_equation(label.y = 1)+labs(x="bsa")
    
    
    # ... Code for final_table ...
    final_table <- bradford_raw[, 4:6]
    names(final_table) <- c("sample", "abs1", "abs2")
    final_table$absorbance <- rowMeans(final_table[, c("abs1", "abs2")], na.rm = TRUE)
    final_table <- final_table[, c("sample", "absorbance")]
    #`predicting concentration`
    final_table$concentrations<-predict.lm(brd_line, final_table)  # Replace with your code
    
    # Apply the dilution
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
    names(final_table)<-c("Sample","Absorbance","µg/µL protein","Laemli 4X","RIPA","Volume of Sample","Final Volume")
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
        footnote=paste(as.character(total_laemli), "µL of Laemli 4X will be needed", sep=" "),
        locations = NULL,
        placement = "auto"
      )
    
    output$brd_plot <- renderPlot({
      brd_plot
    })
    
    output$table <- render_gt({
      tbl
    })
  })
  
  # Download the table as PDF
  output$downloadTable <- downloadHandler(
    filename = function() {
      paste("table", Sys.Date(), ".pdf", sep = "_")
    },
    content = function(file) {
      gtsave(tbl, file = file)
    }
  )
}

# Run the application
shinyApp(ui = ui, server = server)
