
############################
## Ch04 data manipulation ##
############################

library(tidyverse)


## Import csv files

airway <- read_csv("./data/ch04_airway.csv")

class(airway)
dim(airway)
colnames(airway)
head(airway)
tail(airway)

str(airway)

View(airway)

# This training data was generated based on :
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1275876

## Select columns

dplyr::select(airway, Ensembl, value, dex, albut, cell, gene, description)
dplyr::select(airway, -avgLength, -description)


## Select columns based on content = Filtering

filter(airway, dex=="trt")
filter(airway, dex=="trt"&albut=="untrt")

## Pipe operator

# What if you wanted to do all the steps at the same time and store it
# in a new variable?
# You can use intermediate steps, nest the function ... or pipes
# Remember

a1 <- 
  airway %>% # pipe operator shortcut is ctrl + shift + m
  filter(dex=="trt"&albut=="untrt") %>% 
  select(airway, Ensembl, value, dex, albut, cell, gene, description)


## Exercise

# using pipes, subset the airway data with the genes that have a value above 0 counts
# and keep only the sample, ensembl, gene and description.






## Add new columns based on existing columns using mutate

# Frequently you’ll want to create new columns based on the values in 
# existing columns, for example to do unit conversions, or to find the 
# ratio of values in two columns. For this we’ll use mutate().

airway %>% 
  mutate(logval = log10(value)) %>% 
  dplyr::select(SampleName, gene, description, logval, value)

table(is.na(airway))
airway_nona <- airway %>% drop_na()

airway_nona %>% 
  mutate(logval = log10(value)) %>% 
  dplyr::select(SampleName, gene, description, logval, value)

airway %>% 
  filter(!is.na(value)) %>% 
  mutate(logval = log10(value)) %>% 
  dplyr::select(SampleName, gene, description, logval, value)

## Summarize, group_by

# Group_by does not perform any data processing, it groups data into subsets
# for example dex will be divided into treated, untreated and NA

# Group_by is often used with summarize

airway %>% 
  filter(!is.na(value)) %>% 
  group_by(dex) %>% 
  summarize(mean_value = mean(value))

# Here we can answer the question : what's the average gene count per treatment


# But we can also ask per species + per sex

airway %>% 
  filter(!is.na(value)) %>% 
  group_by(dex, cell) %>% 
  summarize(mean_value = mean(value))

airway %>% 
  filter(!is.na(value)) %>% 
  group_by(dex, cell) %>% 
  summarize(mean_value = mean(value)) %>% 
  arrange(desc(mean_value))

## Get summary of columns that do not comprise numerical information

airway %>% 
  count(cell, sort = TRUE)

airway %>% 
  count(cell, Experiment)


# Practice for what we learned -------------------------------------------------

# How many genes have been counted more than 100 times across all samples?

airway %>% filter(value>=100) %>% dim()
airway %>% filter(value>=100)

# (the same gene might have been counted more than 100 times in two different
# samples) We can check if it is the case by doing:

temp <- airway %>% filter(value>=100)
length(unique(temp$gene))
table(duplicated(temp$gene))

# How many genes have been counted more than 100 times in the SRS508568 sample?

airway %>% filter(value>=100, Sample == "SRS508568") %>% dim()

# Select the top 15 highest count values in all samples and make a 
# new table with their description, the run they belong to and the gene names

a2 <- 
  airway %>% 
  arrange(desc(value)) %>% 
  head(15) %>% 
  select(Run,description, gene)

a2

# How many different genes are in this table?

table(a2$gene)

# or 

a2 %>% 
  count(gene)

# What's the average gene count per dex treatment
# For the most counted gene and what's its description?

airway %>% 
  filter(!is.na(value), gene=="CEMIP") %>% 
  group_by(dex) %>% 
  summarize(mean_value = mean(value))

# For that gene, do average counts differ per cells and treatment?
# In which cell and treatment (combined) has it been counted the most?

airway %>% 
  filter(!is.na(value), gene=="CEMIP") %>% 
  group_by(dex, cell) %>% 
  summarize(mean_value = mean(value))



