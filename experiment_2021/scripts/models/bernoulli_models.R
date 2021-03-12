#individual level infection risk models using final infection status
#bernoulli--needed if you're going to include spatial information unique to each individual plant...only ~36000 plants.

#LOAD SOURCE CODE---------------------------------

rm(list=ls())
library(tidyverse)
library(cowplot)
theme_set(theme_classic() + theme(plot.title = element_text(hjust = .5))) 

#IMPORT DATA AND SPATIAL DATAFRAME----------------
plant_level_data <- read_csv('experiment_2021/output/plant_level_data.csv')
tray_level_data <- read_csv('experiment_2021/output/tray_level_data.csv')