options(java.parameters = "-Xss2048k")
library(ggpubr)
library(readxl)
library(caret)
library(dplyr)
library(tidyr)
library(plotly)
library(gt)
library(webshot2)
library(shinythemes)

ui <- fluidPage(
  titlePanel("Bradford Analysis for Western Blots"),
  theme = shinytheme("flatly"),  # Change the overall theme to 'flatly'
  navbarPage(
    "WB Volumes Chart",
    theme = shinytheme("flatly"),  # Change the theme for the navbar
    tabPanel("Analysis",
             sidebarLayout(
               sidebarPanel(
                 actionButton("openLink", "Sample.xlsx", class = "btn btn-info"),  # Set button color to info (blue)
                 fileInput("file", "Choose Excel File", accept = c(".xlsx")),
                 sliderInput("micrograms", "Micrograms Sample", min = 0, max = 100, value = 60),
                 textInput("dilution", "Dilution", value = "1.5"),
                 downloadButton("downloadTable", "Download Table", class = "btn btn-success")  # Set button color to success (green)
               ),
               mainPanel(
                 plotOutput("brd_plot"),
                 gt_output("table")
               )
             )
    ),
    tabPanel("Instructions",
             h3("How to Use the App"),
             p("1. Click on the 'Sample.xlsx' button to open a sample excel file as a template."),
             p("2. Use the 'Choose Excel File' button to upload your own Excel file. The data must be located as in Sample.xlsx."),
             p("3. Adjust the 'Micrograms Sample' slider to set the desired value of protein micrograms per well."),
             p("4. Enter the 'Dilution' value to specify the dilution factor used in samples. (e.g If you added 50 µL of Sample into 50 µL of Bradford the dilution will equal 2"),
             p("5. Click on the 'Run Analysis' button to perform the analysis and generate the plot and table."),
             p("6. The plot will show the quality of the Bradford prediction."),
             p("7. The table will display the calculated volumes to load in the running wells."),
             p("8. Use the 'Download Table' button to download the table as a PNG image.")
    ),
    tabPanel("About",
             h3("Bradford Analysis App"),
             p("This app is designed to perform Bradford analysis on protein samples based on spectrophotometric measurements."),
             p("Created by Andrés Gordo for the ILDH team."),
             h2("Copyright"),
             p("© 2023 Andrés Gordo Ortiz. Attribution 4.0 International (CC BY 4.0)")
              
    )
  )
)

# Define the server
server <- function(input, output) {
  
  observeEvent(input$openLink, {
    browseURL("https://view.officeapps.live.com/op/view.aspx?src=https%3A%2F%2Fraw.githubusercontent.com%2FAndresBioMed%2FUCL-ILDH-Results%2Fmain%2FBradford_calc_app%2Fbradfordtest.xlsx&wdOrigin=BROWSELINK")  # Replace with your desired web link
  })
  
  # Reactive value for micrograms_sample
  micrograms_sample <- reactive({
    input$micrograms
  })
  
  # Read the Excel file and generate brd_plot and table
  observeEvent(input$file, {
    bradford_raw <- read_excel(input$file$datapath)
    bradford_train <- data.frame(absorbance = unlist(bradford_raw[1:9, 2:3]), concentration = 0)
    bradford_train[c(2,11), 2] <- 5
    bradford_train[c(3,12), 2] <- 10
    bradford_train[c(4,13), 2] <- 15
    bradford_train[c(5,14), 2] <- 20
    bradford_train[c(6,15), 2] <- 25
    bradford_train[c(7,16), 2] <- 30
    bradford_train[c(8,17), 2] <- 35
    bradford_train[c(9,18), 2] <- 40
    brd_line <- lm(concentration ~ absorbance, bradford_train)
    
    # Create the plot
    brd_plot <- ggplot(bradford_train, aes(x = concentration, y = absorbance)) +
      geom_smooth(method = "lm", formula = y ~ x, se = FALSE) +
      geom_point() +
      stat_cor(label.y = 0.9, label.x = 0.1) +
      stat_regline_equation(label.y = 0.95, label.x = 0.1) +
      labs(x = "Concentration", y = "Absorbance", title = "Bradford Training Data") +
      theme_bw(base_size = 12, base_family = "Arial") +
      theme(
        plot.title = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        legend.position = "bottom"
      )
    
    output$brd_plot <- renderPlot({
      brd_plot
    })
    
    # Generate the final table
    final_table <- bradford_raw[, 4:6]
    names(final_table) <- c("Sample", "abs1", "abs2")
    final_table$absorbance <- rowMeans(final_table[, c("abs1", "abs2")], na.rm = TRUE)
    final_table <- final_table[, c("Sample", "absorbance")]
    final_table$concentrations <- predict(brd_line, final_table)
    
    # Apply dilution and calculate volumes
    dilution <- as.numeric(input$dilution)
    final_table$concentrations <- final_table$concentrations * dilution
    final_table$volume_sample <- micrograms_sample() / final_table$concentrations
    final_table$laemli4x <- 0.25 * max(final_table$volume_sample) / 0.75
    final_table$final_volume <- 0.25 * max(final_table$volume_sample) / 0.75 + max(final_table$volume_sample)
    final_table$ripa <- final_table$final_volume - (final_table$volume_sample + final_table$laemli4x)
    
    # Rearrange the columns and rename them
    final_table <- final_table %>%
      mutate_if(is.numeric, round, digits = 3) %>%
      relocate(Sample, absorbance, concentrations, laemli4x, ripa, volume_sample, final_volume) %>%
      rename(`µg/µL protein` = concentrations, `Laemli Solution` = laemli4x,
             RIPA = ripa, `Volume of Sample` = volume_sample, `Final Volume` = final_volume)
    
    total_laemli <- sum(final_table$`Laemli Solution`)
    
    # Create the gt table
    tbl <- gt(final_table) %>%
      tab_header(
        title = md("**Volumes Chart for Western Blots**"),
        subtitle = md("Everything is in **µL**")
      ) %>%
      tab_style(
        style = cell_text(weight = "bold"),
        locations = cells_body(columns = 1)
      ) %>%
      tab_style(
        style = list(cell_fill(color = "#add8e6"), cell_text(weight = "bold")),
        locations = cells_body(columns = 4:6)
      ) %>%
      tab_style(
        style = list(cell_fill(color = "#ffefb5"), cell_text(weight = "bold")),
        locations = cells_body(columns = 7)
      ) %>%
      tab_footnote(
        footnote = paste(
          round((total_laemli) / 0.75),
          "µL of Laemli Solution will be needed. Add",
          round((total_laemli) / 0.75) * 0.1,
          "µL of Mercapto 10% into",
          round((total_laemli) / 0.75) * 0.9,
          "µL of Laemli 4X",
          sep = " "
        ),
        locations = NULL,
        placement = "auto"
      )
    
    output$table <- render_gt({
      tbl
    })
  })
  
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
      webshot2::webshot(url = tmp_file, file = file, delay = 0)
      
      # Delete the temporary HTML file
      file.remove(tmp_file)
    }
  )
}

# Run the application
shinyApp(ui = ui, server = server)
