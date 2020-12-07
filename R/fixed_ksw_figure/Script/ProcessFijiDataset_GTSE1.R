# Script to process multiple csv files containing output from Fiji
# Written by James Shelford, generalised by Stephen Royle
## Modified for the Fixed GTSE1 KSW expt. The output from Fiji is already bg corrected.

# wd should be fixed_ksw_figure
# setup preferred directory structure in wd
ifelse(!dir.exists("Data"), dir.create("Data"), "Folder exists already")
ifelse(!dir.exists("Output"), dir.create("Output"), "Folder exists already")
ifelse(!dir.exists("Output/Data"), dir.create("Output/Data"), "Folder exists already")
ifelse(!dir.exists("Output/Plots"), dir.create("Output/Plots"), "Folder exists already")
ifelse(!dir.exists("Script"), dir.create("Script"), "Folder exists already")

# select directory that contains the csv files (name of this folder is reused later)
datadir <- rstudioapi::selectDirectory()
my_files <- list.files(datadir, pattern = "*.csv", full.names = TRUE)
my_files_names <- list.files(datadir, pattern = "*.csv")

# Extract the experiment number (will be used later). Useful for combining data from multiple experiments
Experiment_number <- basename(datadir)

# base R is used to pull in all the data
my_matrix <- matrix(0, length(my_files), 18)

# function definition
build_matrix <- function(my_matrix, my_filename, row_number){
  
  # import data
  my_raw_data <- read.csv(file = my_filename, header = TRUE, stringsAsFactors = FALSE)
  
  # take mean column and transpose
  my_data <- subset(my_raw_data, select = Mean)
  my_data <- t(my_data)
  my_matrix[row_number, 1:18] <- my_data[1, 1:18]
  return(my_matrix)
}

# call the function for each file in the list
for(i in 1:length(my_files)){
  my_filename <- my_files[i]
  my_matrix <- build_matrix(my_matrix, my_filename, i)
}

# Generating the mean values and ratio
reference_spindle_matrix <- matrix(0, length(my_files), 3)
reference_spindle_matrix[1:length(my_files), 1:3] <- my_matrix[1:length(my_files), c(1:3)]

GFP_spindle_matrix <- matrix(0, length(my_files), 3)
GFP_spindle_matrix[1:length(my_files), 1:3] <- my_matrix[1:length(my_files), c(4:6)]

POI_spindle_matrix <- matrix(0, length(my_files), 3)
POI_spindle_matrix[1:length(my_files), 1:3] <- my_matrix[1:length(my_files), c(7:9)]

reference_cytoplasm_matrix <- matrix(0, length(my_files), 3)
reference_cytoplasm_matrix[1:length(my_files), 1:3] <- my_matrix[1:length(my_files), c(10:12)]

GFP_cytoplasm_matrix <- matrix(0, length(my_files), 3)
GFP_cytoplasm_matrix[1:length(my_files), 1:3] <- my_matrix[1:length(my_files), c(13:15)]

POI_cytoplasm_matrix <- matrix(0, length(my_files), 3)
POI_cytoplasm_matrix[1:length(my_files), 1:3] <- my_matrix[1:length(my_files), c(16:18)]

reference_spindle_means <- rowMeans(reference_spindle_matrix, na.rm = TRUE)
reference_cytoplasm_means <- rowMeans(reference_cytoplasm_matrix, na.rm = TRUE)
POI_spindle_means <- rowMeans(POI_spindle_matrix, na.rm = TRUE)
POI_cytoplasm_means <- rowMeans(POI_cytoplasm_matrix, na.rm = TRUE)
GFP_spindle_means <- rowMeans(GFP_spindle_matrix, na.rm = TRUE)
GFP_cytoplasm_means <- rowMeans(GFP_cytoplasm_matrix, na.rm = TRUE)

# ratio calculation
reference_spindle_ratio <- reference_spindle_means / reference_cytoplasm_means
POI_spindle_ratio <- POI_spindle_means / POI_cytoplasm_means
GFP_spindle_ratio <- GFP_spindle_means / GFP_cytoplasm_means

# Adding these values to the matrix
my_matrix <- cbind(my_matrix, reference_spindle_means, POI_spindle_means, GFP_spindle_means, reference_cytoplasm_means, POI_cytoplasm_means, GFP_cytoplasm_means, reference_spindle_ratio, POI_spindle_ratio, GFP_spindle_ratio)

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

# load external look-up table
look_up_table <- read.table("Data/lookup.csv", header = TRUE, stringsAsFactors = FALSE, sep = ",")

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

# Add the experiment number to the df
df1$Experiment_number <- Experiment_number

cellLineName <- basename(dirname(datadir))
pathToDataSave <-  paste0("Output/Data/", cellLineName, "/")
ifelse(!dir.exists(pathToDataSave), dir.create(pathToDataSave), "Folder exists already")

# Now save the dataset the name will be the folderName that the csvs were in
folderName <- basename(datadir)
fileName <- paste0(pathToDataSave, "dataframe_", folderName, ".rds")
saveRDS(df1, file=fileName)

