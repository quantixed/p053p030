## Fixed knocksideways data visualisation
## Load in dataframes produced from quantifying individual experiments, then plot graph with stats.
## Original script by James Shelford
## Add-ons and generalisation by Stephen Royle

require(ggplot2)
library(dplyr)
library(multcomp)
library(scales)
library(cowplot)

# make directory for output if it doesn't exist
if (dir.exists("output")==FALSE) dir.create("output")

# Load the dataframes for each experiment
datadir <- rstudioapi::selectDirectory()
my_files <- list.files(datadir,pattern='*.rds',full.names = TRUE)
combined_df <- bind_rows(lapply(my_files, readRDS))
combined_df$Category <- as.factor(combined_df$Category)

# How many cells in each condition?
summary(combined_df$Category)

# To plot in the correct order, we first need the look-up table
look_up_table <- read.table("lookup.csv", header = TRUE, stringsAsFactors = F, sep = ",")
combined_df$Category <- factor(combined_df$Category, levels = look_up_table$Search_category )

# Subsetting the data to remove tubulin and PI3KC2A
theNames <- rev(names(combined_df))
theNames <- theNames[1:4]
combined_df <- subset(combined_df, Category != "tubulin control", select = theNames)
combined_df <- subset(combined_df, Category != "PI3KC2A control", select = theNames)
combined_df <- subset(combined_df, Category != "tubulin rapamycin", select = theNames)
combined_df <- subset(combined_df, Category != "PI3KC2A rapamycin", select = theNames)
names(combined_df) <- c("Category","blind_list","GFP_spindle_ratio","POI_spindle_ratio")

## Statistics
the_ANOVA <- aov(GFP_spindle_ratio ~ Category, combined_df)
summary(the_ANOVA)
plot(the_ANOVA,1)
plot(the_ANOVA,2)
summary(glht(the_ANOVA, linfct = mcp(Category='Tukey')))

## Generating the plot
# we need to know the name of the cell line
theCellLine <- basename(datadir)

# find the min and max of the dataframe
loVal <- min(combined_df[,3:4], na.rm=T)
hiVal <- max(combined_df[,3:4], na.rm=T)

# we want symmetry about 1
findAxisLimit <- function(val1,val2){
  highVal <- max(1/val1,val2)
  highVal <- ceiling(log2(highVal))
  highVal <- 2^highVal
  return(highVal)
}
axVal <- findAxisLimit(loVal,hiVal)

# function to generate the plots
makeTheScatterPlot <- function(preData, postData, yLab, xLab) {
  ggplot(filter(combined_df, Category == preData | Category == postData),
         aes(x=GFP_spindle_ratio, y=POI_spindle_ratio, color=Category, alpha=0.5)) +
    geom_point() +
    scale_x_continuous(trans='log2', limits = c(1/axVal,axVal), breaks = trans_breaks("log2", function(x) 2^x), labels = trans_format("log2", math_format(.x))) +
    scale_y_continuous(trans='log2', limits = c(1/axVal,axVal), breaks = trans_breaks("log2", function(x) 2^x), labels = trans_format("log2", math_format(.x))) +
    labs(y = yLab, x = xLab) +
    theme(legend.position = 'none')
}

# make each of the four plots
p1 <- makeTheScatterPlot("CHC control","CHC rapamycin", "CHC", theCellLine)
p2 <- makeTheScatterPlot("TACC3 control","TACC3 rapamycin", "TACC3", theCellLine)
p3 <- makeTheScatterPlot("chTOG control","chTOG rapamycin", "chTOG", theCellLine)
p4 <- makeTheScatterPlot("GTSE1 control","GTSE1 rapamycin", "GTSE1", theCellLine)

# arrange the plots, display and save as PDF
all_scatter_plots <- plot_grid(p1,p2,p3,p4, rel_widths = c(1, 1), rel_heights = c(1,1)) + theme(aspect.ratio=1)
all_scatter_plots
ggsave("./output/all_scatter_plots.pdf", plot = all_scatter_plots, width = 120, height = 120, units = 'mm', useDingbats = FALSE)
#ggsave("./output/all_scatter_plots.png", plot = all_scatter_plots, dpi = 300)
