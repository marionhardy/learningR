
############################
## Ch03 data manipulation ##
############################

# Correction of previous exercises

## Import csv files

surveys <- read.csv("./data/ch03_surveys_data.csv")

class(surveys)
dim(surveys)
colnames(surveys)

head(surveys)
tail(surveys)

str(surveys)

View(surveys)

## Using brackets to subset dataframes

# Removing the first column

s1 <- surveys[,-1]

head(s1) # first column has been removed
s1[c(1:5),] # the equivalent of the line above

# Subset the dataframe in other ways

s1[,c("record_id", "species_id","sex")]
s1[,-c("plot_type")] # will not work
s1[,-13] # will work
s1[,-ncol(surveys)]
s1[,colnames(surveys)!="plot_type"]

s_1977 <- s1[s1$species==1977,]
table(s1$year)
table(s_albigula$year)

s_bigfoot <- s1[s1$hindfoot_length>= 28,]

s_F_al <- s1[s1$sex=="F"&s1$species=="albigula",]

table(s_F_al$sex)
table(s_F_al$species)

FALSE&TRUE # it's and
TRUE&TRUE
FALSE&FALSE

s_1977_78 <- s1[s1$year==1977|s1$year==1978,]

table(s_1977_78$year)
 

TRUE|FALSE # it's a and/or
FALSE|TRUE
TRUE|TRUE
FALSE|FALSE

## Save your data

write.csv(s_F_al,'./data/surveys_female_and_albigula.csv')


## Now susbet surveys and select the mice that are of the albigula,
# flavus and savannarum species with a hindfoot length of above 15cm
# remove the first column and make a new column with the hindfoot
# length of the mice converted into millimiters

s2 <- surveys[s1$species=="savannarum"|
                s1$species=="albigula"|
                s1$species=="flavus",-1]

s2 <- s2[s2$hindfoot_length>=25,] # annoying hindfoot data missing

s2$hindfoot_mm <- s2$hindfoot_length*10

head(s2)

## Install tidyverse package

# install.packages("tidyverse")
library(tidyverse)

# This tool will allow us to make operations on dataframe in a way more
# efficient and readable way
# the tidyverse packages contains (among other things)
# the dplyr, ridyr, ggplot2 packages

# dplyr is a package for making tabular data (dataframes) much
# easier to manipulate. It contains functions that are verbs
# which make data manipulation much more grammaticaly readable

# tidyr is a package that helps you easily convert data formats
# it makes data tidier and helps you handle missing values
# drop_na, fill, replace_na ...

# ggplot2 is a package used to make plots and visualize your data 
# with an intuitive grammar to do so and multiple plotting functions
# violin plots, density plots and much much more..

# Let's start again from our surveys dataframe

surveys <- read_csv('./data/ch03_surveys_data.csv')
str(surveys)

## Select columns

select(surveys, record_id, year, sex, species)

select(surveys, -genus, -plot_type)

## Select columns based on content = Filtering

filter(surveys, sex=="M")
filter(surveys, sex=="M"& year == 1977)

## Pipe operator

# What is you wanted to do all the steps at the same time?
# You can use intermediate steps, nest the function ... or pipes
# Remember

s2 <- surveys[s1$species=="savannarum"|
                s1$species=="albigula"|
                s1$species=="flavus",-1]

# Our tibble does not have the first column duplicated

s2 <- s2[s2$hindfoot_length>=25,] 

# Now you can select the same things but by doing:


s3 <- select(filter(surveys,hindfoot_length>=25), species, hindfoot_length, sex)



s3 <- 
surveys %>% # pipe operator shortcut is ctrl + shift + m
  filter(hindfoot_length>=25, species %in% c("savannarum","albigula","flavus")) %>% 
  select(species, sex, hindfoot_length)


## Exercise

# using pipes, subset surveys with the data form the years 1977 and 1978
# and retain only the columns species_id, year, sex, hindfoot_length






## Add new columns using mutate

# Frequently you’ll want to create new columns based on the values in 
# existing columns, for example to do unit conversions, or to find the 
# ratio of values in two columns. For this we’ll use mutate().

surveys %>% 
  mutate(hindfoot_mm = hindfoot_length*10) %>% 
  select(species_id, sex, hindfoot_mm, hindfoot_length)


surveys_nona <- surveys %>% drop_na()

surveys_nona %>% 
  mutate(hindfoot_mm = hindfoot_length*10) %>% 
  select(species_id, sex, hindfoot_mm, hindfoot_length)


surveys %>% 
  filter(!is.na(hindfoot_length)) %>% 
  mutate(hindfoot_mm = hindfoot_length*10) %>% 
  select(species_id, sex, hindfoot_mm, hindfoot_length)


surveys %>% 
  filter(!is.na(hindfoot_length)) %>% 
  mutate(hindfoot_mm = hindfoot_length*10) %>% 
  select(species_id, sex, hindfoot_mm, hindfoot_length) %>% 
  head()



## Exercise: Create a new dataframe that contains only the sex
# species id, species, year and hindfoot length columns. Add a new
# column that contains the log() of the hindfoot length. This dataframe
# must only contain data from male mice who must be associated to an
# existing record_id






## Summarize, group_by

# Group_by does not perform any data processing, it groups data into subsets
# for example sex will be subset into male, female, NA
# operations can be carried on, on those temporary grouped subsets

# Group_by is often used with summarize


surveys %>% 
  group_by(species) %>% 
  summarize(mean_hindfoot = mean(hindfoot_length))

# Here we can answer the question : what's the average hindfoot
# length per species

# But we can also ask per species + per sex

surveys %>% 
  group_by(species, sex) %>% 
  filter(!is.na(hindfoot_length)) %>% 
  summarize(mean_hindfoot = mean(hindfoot_length))

surveys %>% 
  group_by(species, sex) %>% 
  filter(!is.na(hindfoot_length)) %>% 
  summarize(mean_hindfoot = mean(hindfoot_length),
            median_hindfoot = median(hindfoot_length)) %>% 
  arrange(mean_hindfoot)

surveys %>% 
  group_by(species, sex) %>% 
  filter(!is.na(hindfoot_length)) %>% 
  summarize(mean_hindfoot = mean(hindfoot_length),
            median_hindfoot = median(hindfoot_length)) %>% 
  arrange(desc(mean_hindfoot))


surveys %>% 
  group_by(species) %>% 
  count(n())

surveys %>% 
  count(species, sort = TRUE)

surveys %>% 
  count(species, year)






