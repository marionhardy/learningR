#########################
## Data visualization ##
########################

library(tidyverse)
# tidyverse also gives access to premade datasets such as iris, mtcars etc

# BiocManager::install("ggplot2")
library(ggplot2)
# the package commonly used by r programmers to plot figures


## Load the data we will use for this course

data("airquality") # double click on the <Promise> which appeared in the 
                  # global environment

# This dataset records the airquality in New York from May 1, 1973 to September 30, 1973.

table(is.na(airquality))
colSums(is.na(airquality))

# If we want to plot this information, we can use ggplot()

ggplot(data = airquality) # a plot appears in the plot window

ggplot(data = airquality,
       aes(x = Month, y = Temp)) # in the aesthetic argument, you can input
                                # what it is you want to plot on the axes

ggplot(data = airquality,
       aes(x = Month, y = Temp))+ # your plot now contains dots representing the data
  geom_point()

ggplot(data = airquality,
       aes(x = Month, y = Temp, color = Day))+ # you can color by a variable
  geom_point()
# Two things seems weird though, 
# 1. there do not seem to be 31 dots for each month
# 2. the colors are not separate, it's a gradient
# what happened?

#1. plot the same points but making them translucent to a certain degree
ggplot(data = airquality,
       aes(x = Month, y = Temp))+
  geom_point(alpha = 0.1)

#2. Check the class of Day -> which is it?
class(airquality$Day) # integer is similar to numeric
# what would you change it into to make R understand that they are distinct from 
# each other?

airquality$Day <- as.character(airquality$Day)

ggplot(data = airquality,
       aes(x = Month, y = Temp, color = Day))+
  geom_point() # colors are now distinct but what about overlapping data?

# One solution is to keep the alpha = 0.1 parameter
# One solution is to slightly offset each point form another using jitter

ggplot(data = airquality,
       aes(x = Month, y = Temp, color = Day))+
  geom_point()+
  geom_jitter() # However, Distinction between months is less clean

# Another solution would be to separate the plot per months

ggplot(data = airquality,
       aes(x = Day, y = Temp, color = Day))+
  geom_point()+
  facet_grid(.~ Month) 

ggplot(data = airquality,
       aes(x = Day, y = Temp, color = Day))+
  geom_point()+
  facet_grid(Month~.) # This facetting makes the graph more readable

ggplot(data = airquality,
       aes(x = Day, y = Temp, color = Ozone, size = Wind))+
  geom_point()+
  facet_grid(Month~.) # you can add lots of information in this graph

# You can use other representation than points
# and you can combine ggplot with the selection from your %>% 
# for example, plot the temperatures from only the month of May using lines

airquality %>% 
  filter(Month==5) %>% 
  ggplot(aes(x=Day, y=Temp))+
  geom_line()

# Now plot every month

airquality %>% 
  ggplot(aes(Day, Temp, group = Month))+
  geom_line()+
  facet_grid(Month~.)

# Something is wrong, can you see what?

# By changing days into characters, we loose their numerical increasing order
# a solution to this which would also make us keep the distinct colors and
# not reverse to a gradient for coloring -> transform Days into factors

airquality$Day <- as.factor(as.numeric(airquality$Day)) # which have levels

airquality %>% 
  ggplot(aes(Day), Temp, group = Month)+
  geom_line()+
  facet_grid(Month~.)

# Everything is customizable

airquality %>% 
  ggplot(aes(Day, Temp, group = Month))+
  geom_line(color = "blue", size = 3, alpha = 0.5)+
  facet_grid(Month~.)+
  theme_bw()

# To come back to our representations of temperature per days per months,
# a truly more helpful representation of the distribution is boxplots
# and violinplots which are accessible in ggplot

ggplot(data = airquality,
       aes(x = Month, y = Temp, group = Month))+
  geom_boxplot()+
  theme_bw()

ggplot(data = airquality,
       aes(x = Month, y = Temp, group = Month, color = Month))+
  geom_violin()+
  theme_bw()

# We changed Day but did not change Month! Be careful to inspect your data 
# thoroughly

ggplot(data = airquality,
       aes(x = Month, y = Temp, group = Month, color = Month))+
  geom_violin()+
  theme_classic()+
  labs(title="Temperature per month", subtitle = "221019")

ggplot(data = airquality,
       aes(x = Month, y = Temp, group = Month, color = Month))+
  geom_violin(trim = F)+
  geom_jitter(position = position_jitter(0.2), alpha = 0.5)+
  geom_boxplot(width = 0.1)+
  geom_hline(yintercept = 65, linetype = 4)+
  theme_classic()+
  labs(title="Temperature per month", subtitle = "221019")


## Exercise

  # Plot temperature and Solar rays per month. Does it seem like
  # Temperature and solar radiation are correlated?

airquality %>% 
  ggplot(aes(Temp, Solar.R))+
  geom_point()+
  facet_grid(Month~.)

  # Change months into their given names and add, in your chosen
  # representation, a way to also see wind strength
  # add a title to your graph

airquality <- airquality %>% 
  arrange(Month)

airquality$Months <- as.factor(c(rep('May', 31), rep("June", 30), rep("July", 31),
                          rep("August", 31), rep("September", 30)))

airquality$Months <- ordered(airquality$Months, levels = c("May","June","July","August","September"))

airquality %>% 
  ggplot(aes(Temp, Solar.R, size = Wind))+
  geom_point(color = "dark blue")+
  facet_grid(Months~.)+
  theme_light()+
  labs(title = "Temperature and solar radiation per months")

## To export your plots, you can use ggsave()

ggsave(filename = "./figures/solar_per_temp", plot = last_plot(), device = png, dpi = 400)

# Now run the next line of code:

airquality$Measurer <- c(rep("Laudna",55),rep("Imogen", 65), rep("Ashton", 33))

# Check if the three people collecting samples usually have similar measure values

airquality %>% 
  ggplot(aes(Measurer, Temp, color = Measurer))+
  geom_boxplot()+
  theme_light() # Looks like Laudna has more variation in her measures
              # and they seem to be overall lower

# Check if this is a trend each months

airquality %>% 
  ggplot(aes(Measurer, Temp, color = Measurer))+
  geom_violin()+
  facet_grid(.~Months)
  theme_light()
  
# Now who is the problem that is highlighted by your graph?
  
airquality %>% 
  ggplot(aes(Temp, Solar.R, size = Wind))+
  geom_point(color = "dark blue")+
  facet_grid(Months~Measurer)+
  theme_light()+

# How would you verify whether this difference is due to anything else?

airquality %>% 
  ggplot(aes(Temp, Solar.R, size = Wind, color = Ozone))+
  geom_point()+
  facet_grid(Months~Measurer)+
  theme_light()

## RNAseq data for volcano plot + highlight per filtering

rna <- read_csv("./data/res_tbl.csv")
str(rna)

# To plot RNAseq data, a violinplot is very common

rna %>%
  filter(!is.na(padj)) %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(size = 0.5) 

rna %>%
  filter(!is.na(padj)) %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(size = 0.5)+
  geom_hline(yintercept = -log10(0.05)) +
  geom_vline(xintercept = 1) +
  geom_vline(xintercept = -1)+
  labs(title = "Selected cells T96 vs T48")+
  theme_bw()

rna %>%
  filter(!is.na(padj)) %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj),
              color = padj < 0.05 & abs(log2FoldChange) > 1)) +
  scale_colour_manual(values = c("gray", "firebrick3")) +
  geom_point(size = 0.5) +
  geom_hline(yintercept = -log10(0.05)) +
  geom_vline(xintercept = 1) +
  geom_vline(xintercept = -1)+
  labs(title = "Selected cells T96 vs T48")+
  theme_bw()

rna %>%
  filter(!is.na(padj)) %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj),
             color = padj < 0.05 & abs(log2FoldChange) > 1, 
             label = ifelse(padj<0.05&log2FoldChange>=1|padj<0.05&log2FoldChange<=-1,
                            as.character(gene),''))) +
  scale_colour_manual(values = c("gray", "firebrick3")) +
  geom_point(size = 0.5) +
  geom_hline(yintercept = -log10(0.05)) +
  geom_vline(xintercept = 1) +
  geom_vline(xintercept = -1)+
  geom_text_repel(max.overlaps = 10)+
  labs(title = "Selected cells T96 vs T48")+
  theme_bw()


## For your information, base R also has plot functions which function
# very well and can be used when doing quality control where there is no
# need for fancy colors and titles


hist(rna$padj)

plot(rna$log2FoldChange, -log10(rna$padj))
abline(a = -log10(0.05), b= 0, v = c(1,-1))

boxplot(airquality$Temp ~ airquality$Measurer)

dev.null()




