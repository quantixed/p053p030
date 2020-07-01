# Script to process multiple csv files containing output from Fiji
# Written by James Shelford, generalised by Stephen Royle

# log.txt needs to be in the "data" folder

# Select folder where analysis takes place
analysisdir <- rstudioapi::selectDirectory()
setwd(analysisdir)

# make directory for output if it doesn't exist
if (dir.exists("output") == FALSE) dir.create("output")

# select directory that contains the csv files (name of this folder is reused later)
datadir <- rstudioapi::selectDirectory()
my_files <- list.files(datadir, pattern = "*.csv", full.names = TRUE)
my_files_names <- list.files(datadir, pattern = "*.csv")

# base R is used to pull in all the data
my_matrix <- matrix(0, length(my_files), 12)

# function definition
build_matrix <- function(my_matrix, my_filename, row_number){
  
  # import data
  my_raw_data <- read.csv(file = my_filename, header = TRUE, stringsAsFactors = FALSE)
  
  # take mean column and transpose
  my_data <- subset(my_raw_data, select = Mean)
  my_data <- t(my_data)
  my_matrix[row_number, 1:12] <- my_data[1, 1:12]
  return(my_matrix)
}

# call the function for each file in the list
for(i in 1:length(my_files)){
  my_filename <- my_files[i]
  my_matrix <- build_matrix(my_matrix, my_filename, i)
}

# Generating the values and ratio
#In this matrix, the first 3 values are the MT intensity and the next 3 values are the corresponding bg values. Tubulin first and then POI

Tubulin_MT_matrix <- matrix(0,length(my_files),3)
Tubulin_cytosol_matrix <- matrix(0, length(my_files),3)
Tubulin_MT_matrix[1:length(my_files),1:3] <- my_matrix[1:length(my_files), c(1:3)]
Tubulin_cytosol_matrix[1:length(my_files),1:3] <- my_matrix[1:length(my_files), c(7:9)]

POI_MT_matrix <- matrix(0,length(my_files),3)
POI_cytosol_matrix <- matrix(0, length(my_files),3)
POI_MT_matrix[1:length(my_files),1:3] <- my_matrix[1:length(my_files), c(4:6)]
POI_cytosol_matrix[1:length(my_files),1:3] <- my_matrix[1:length(my_files), c(10:12)]

Tubulin_MT_means <- rowMeans(Tubulin_MT_matrix, na.rm = TRUE)
Tubulin_cytosol_means <- rowMeans(Tubulin_cytosol_matrix, na.rm = TRUE)

POI_MT_means <- rowMeans(POI_MT_matrix, na.rm = TRUE)
POI_cytosol_means <- rowMeans(POI_cytosol_matrix, na.rm = TRUE)

Tubulin_MT_ratio <- Tubulin_MT_means / Tubulin_cytosol_means
POI_MT_ratio <- POI_MT_means / POI_cytosol_means


# Adding these values to the matrix
my_matrix <- cbind(my_matrix, Tubulin_MT_means, Tubulin_cytosol_means, POI_MT_means, POI_cytosol_means, Tubulin_MT_ratio, POI_MT_ratio)

# Make list of the names of the blinded filenames with *.csv removed put it all together in data frame
blind_list <- gsub(".csv","", my_files_names)
df1 <- as.data.frame(my_matrix)
df1$blind_list <- blind_list

# load the log.txt file
logfile_path <- paste0(datadir,"/log.txt")
blind_log <- read.table(logfile_path, header = TRUE)

# function to find partial strings in a column and classify them
add_categories = function(x, patterns, replacements = patterns, fill = NA, ...) {
  stopifnot(length(patterns) == length(replacements))
  ans = rep_len(as.character(fill), length(x))    
  empty = seq_along(x)
  for(i in seq_along(patterns)) {
    greps = grepl(patterns[[i]], x[empty], ...)
    ans[empty[greps]] = replacements[[i]]  
    empty = empty[!greps]
  }
  return(ans)
}

#Load the look-up table

look_up_table <- read.table('lookup.csv', header = TRUE, stringsAsFactors = FALSE, sep = ",")

# add a new column to dataframe where categories are defined by searching original name for partial strings
blind_log$Category <- add_categories(blind_log$Original_Name,
                                     look_up_table$Search_name,
                                     look_up_table$Search_category,
                                     "NA", ignore.case = TRUE)
# Now we have a dataframe that can be used to lookup the real values

# This line looks up the correct Category from the blind_log
df1$Category <- with(blind_log,
                     Category[match(df1$blind_list,
                                    Blinded_Name)])
# needs to be done like this because
# a) blind_log is in a random order
# b) your list of *.csv names could be in any order (although they're probably sorted alphanumerically)

# Now save the dataset the name will be the folderName that the csvs were in
folderName <- basename(analysisdir)
fileName <- paste0("./output/dataframe_", folderName, "_interphase.rds")
saveRDS(df1, file=fileName)
