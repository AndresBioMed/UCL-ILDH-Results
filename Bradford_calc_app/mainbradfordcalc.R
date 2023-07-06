library(ggpubr)
library(readxl)
library(caret)
library(dplyr)
library(tidyr)
library(plotly)
library(gt)
#reading the bradford excel
bradford_raw <- read_excel("C:\\Users\\andre\\OneDrive - unizar.es\\Laboratorio\\UCL ILDH 2023\\GitHub Repository UCL\\UCL-ILDH-Results\\bradfordtest.xlsx", col_names = FALSE)
bradford_train<-data.frame(absorbance=unlist(bradford_raw[1:9,2:3]), concentration=0)
bradford_train[c(2,11), 2]<- 5
bradford_train[c(3,12), 2]<- 10
bradford_train[c(4,13), 2]<- 15
bradford_train[c(5,14), 2]<- 20
bradford_train[c(6,15), 2]<- 25
bradford_train[c(7,16), 2]<- 30
bradford_train[c(8,17), 2]<- 35
bradford_train[c(9,18), 2]<- 40
#bradford training
brd_line <- lm(concentration~absorbance, bradford_train)
brd_plot<-ggplot(bradford_train, aes(x=concentration,y=absorbance))
brd_plot+geom_smooth(method = "lm", formula = y ~ x) +geom_point() + stat_cor(label.y = 0.9)+ 
  stat_regline_equation(label.y = 1)+labs(x="bsa")

#data strapping
final_table <- bradford_raw[, 4:6]
names(final_table) <- c("sample", "abs1", "abs2")
final_table$absorbance <- rowMeans(final_table[, c("abs1", "abs2")], na.rm = TRUE)
final_table <- final_table[, c("sample", "absorbance")]
#`predicting concentration`
final_table$concentrations<-predict.lm(brd_line, final_table)
#applying the dilution 
dilution=1.5
micrograms_sample=60
final_table$concentrations<-final_table$concentrations*dilution
#f
final_table$volume_sample<-micrograms_sample/final_table$concentrations
final_table$laemli4x<-0.25*max(final_table$volume_sample)/0.75
final_table$final_volume<-0.25*max(final_table$volume_sample)/0.75+max(final_table$volume_sample)
final_table$ripa<-final_table$final_volume-(final_table$volume_sample+final_table$laemli4x)
final_table<-final_table %>% 
  mutate_if(is.numeric, round, digits=3)%>%
  relocate(sample, absorbance, concentrations, laemli4x, ripa, volume_sample, final_volume)

total_laemli<-sum(final_table$laemli4x)
tbl <- gt(final_table)

# Apply styling options
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
  )# Apply to the last four columns
  
