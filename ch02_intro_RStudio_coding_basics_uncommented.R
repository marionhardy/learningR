
###################################################
## Chapter 2 : The basics of RStudio and coding ##
##################################################

### Good practices

## Create a project: here -> LearningR.proj

getwd()

# Create the generally useful subfolders: data, figures
# Learn to comment your code using the pound (#) sign, it's useful and a must
# for reproducibility

### Start creating objects and using the console or/and a script

5+78
198/3

target_gene <- 6 
gene_all = 8 

print(gene_all)

gene_all

(target_gene/gene_all)*100

### Types of objects 

## The vector : The most common data structure

gene_freq <- c(0.01,20,65,0.3)

gene_names <- c("PRKAA1","MTOR","ULK1","SLC1A4")

gene_names <- c(PRKAA1,MTOR,ULK1,SLC1A4)


######### Data types in R ###########

# Characters: "MTOR", "10", "TRUE"
# Numeric : 2, 15.5, 120.2354
# Logical : TRUE FALSE

#####################################

class(gene_names) 

x <- c("FR", 10, "BE")
x[1] 
x[2]


class(x[2]) 


# Try seeing what will happen in each of those examples

num_chr <- c(1, 2, 3, "a")
num_logi <- c(1, 2, 3, TRUE)
chr_logi <- c("a", "b", "c", TRUE)
tricky <- c(1, 2, 3, "4")


x[c(1,3)]
x[c(3,1)]
x[-2] 
x[-c(2,3)] 

x <- c(x,"This","Is","Longer")

sort(x)


## The list (not Schindler's)

x <- list(1,"PRKAA2", TRUE)
x
class(x[1])
class(x[[1]]) 
class(x[[2]])

## The matrix : Reloaded

m <-matrix(1:6, nrow = 2, ncol = 3)
m
dim(m)
nrow(m)
ncol(m)


x <- c(1,2)
y <- c(3,4)
z <- c(5,6)
m1 <- cbind(x,y,z)
?cbind

# or

x <- c(1,3,5)
y <- c(2,4,6)

m2 <- rbind(x,y)

m==m1
m1==m2


x <- c(1,3,5)
y <- c(2,4,6,7)

m3 <- rbind(x,y)

length(x)
length(y)
length(x)==length(y)


m[1] 
m[1,]
m[,1]

# Try accessing the value in the second row and third column



## The dataframe


set.seed(2)
data <- data.frame(ID = 2355:2361,
                   name = c("PRKAA1", "CFTR", "EGFR","MYC","RAS","P53","CD14"),
                   expr = round(runif(7,0,1), 2))

data

rownames(data)
colnames(data)

data$ID
data$name

sort(data)
data[,data$expr>0.5]

## Additional exercises

# Create a vector of 7 random pvalues between 0.0001 and 0.01 and add it to the "data" dataframe
# hint: create a vector and fill it with the runif() function


# If not already done, call the column you added "pvalue"


# Change the rownames of the "data" dataframe to the IDs of the genes

























