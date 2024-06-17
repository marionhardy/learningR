
###################################################
## Chapter 2 : The basics of RStudio and coding ##
##################################################

### Good practices

## Create a project: here -> LearningR.proj

getwd() # tells you which folder R can see and search from, shoudl be LearningR

# Create the generally useful subfolders: data, figures
# Learn to comment your code using the pound (#) sign, it's useful and a must
# for reproducibility

### Start creating objects and using the console and/or a script

# You can get output from R by simply doing math:

5+78
198/3

# You can also store this data, which will be useful for the future

target_gene <- 6 
gene_all = 8 


print(gene_all)
# or
gene_all

# You can now play with your variable 

(target_gene/gene_all)*100

### Types of objects 

## The vector : The most common data structure

gene_freq <- c(0.01,20,65,0.3)

# A vector can also contain characters

gene_names <- c("PRKAA1","MTOR","ULK1","SLC1A4")
gene_names <- c(PRKAA1,MTOR,ULK1,SLC1A4)


######### Data types in R ###########

# Characters: "MTOR", "10", "TRUE"
# Numeric : 2, 15.5, 120.2354
# Logical : TRUE FALSE

#####################################

class(gene_names) # to inspect the class of elements

# So what happens if you mix them?
x <- c("FR", 10, "BE")
x[1] # to access the first element of the vector use [1]
x[2]

# NB: MATLAB and R start counting from 1
# But C/C++, Python and Java start counting from 0

class(x[2]) # The 10 as been transformed into a "10" which is a character


# Try seeing what will happen in each of those examples

num_chr <- c(1, 2, 3, "a")
num_logi <- c(1, 2, 3, TRUE)
chr_logi <- c("a", "b", "c", TRUE)
tricky <- c(1, 2, 3, "4")

class(num_chr)
class(num_logi)
class(chr_logi)
class(tricky)

# this class coercion has a logic:
# logical → numeric → character ← logical 

# If you wanted to select more elements of the vector, you could use

x[c(1,3)]
x[c(3,1)]
x[-2] 
x[-c(2,3)] 

# If you want to add elements to a vector

x <- c(x,"This","Is","Longer")

# If you want to sort them by alphabetical order

sort(x)

## The list (not Schindler's)

# Behave like containers. Unlike vectors, they can contain
# different data types

x <- list(1,"PRKAA2", TRUE)
x
class(x[1])
class(x[[1]]) 
class(x[[2]])

## The matrix : Reloaded

# Matrices are an extension of vectors. They are vectors with dimensions.
# Meaning that they only tolerate one data type at the same time
# You can check their dimensions by using dim()

m <-matrix(1:6, nrow = 2, ncol = 3)
m
dim(m)
nrow(m)
ncol(m)

# Could also have been created by doing

x <- c(1,2)
y <- c(3,4)
z <- c(5,6)
m1 <- cbind(x,y,z)
?cbind # if you don't know what a function does, this is how you access its documentation

# or

x <- c(1,3,5)
y <- c(2,4,6)

m2 <- rbind(x,y)

# We can check whether or not they are similar by using the ==

m==m1
m1==m2

# If we had tried

x <- c(1,3,5)
y <- c(2,4,6,7)

m3 <- rbind(x,y)

length(x)
length(y)
length(x)==length(y)

# To access the first row of the matrix

m[1]
m[1,]

# To access the first column of the matrix

m[,1]

# Try accessing the value in the second row and third column

m[2,3]


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


# If not already done, 
# call the column you added "pvalue"
# and change the rownames of the "data" dataframe to the IDs of the genes
# remove the ID column



























