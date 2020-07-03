## GTSE1 LIDL mutant data visualisation
## Load in dataframes produced from quantifying individual experiments, then plot graph with stats.
## Original script by James Shelford
## Add-ons and generalisation by Stephen Royle

require(ggplot2)
require(ggbeeswarm)
library(dplyr)
library(multcomp)
library(scales)
library(cowplot)

# make directory for plots if it doesn't exist (it should)
ifelse(!dir.exists("Output"), dir.create("Output"), "Folder exists already")
ifelse(!dir.exists("Output/Plots"), dir.create("Output/Plots"), "Folder exists already")

# Load the dataframes for all experiments for one cell line
# these are rds files found in Output/Data
datadir <- "Output/Data"
my_files <- list.files(datadir,pattern='*.rds',full.names = TRUE)
combined_df <- bind_rows(lapply(my_files, readRDS))
combined_df$Category <- as.factor(combined_df$Category)

# Subset to remove NAs
theNames <- rev(names(combined_df))
theNames <- theNames[1:4]
combined_df <- subset(combined_df, Category != "NA", select = theNames)

# How many cells in each condition?
summary(combined_df$Category)

# To plot in the correct order, we first need the look-up table
look_up_table <- read.table("Data/lookup.csv", header = TRUE, stringsAsFactors = F, sep = ",")
combined_df$Category <- factor(combined_df$Category, levels = look_up_table$Search_category )

## Generating the plots

combined_clathrin_plot <- ggplot(filter (combined_df), aes(x = Category,y = clathrin_spindle_ratio, colour = "#00A651")) +
  scale_colour_manual(values = "#00A651") +
  geom_hline(yintercept =1, linetype='dashed', colour='black') +
  geom_quasirandom(alpha = 0.5, stroke = 0) +
  stat_summary(fun.data = mean_se, geom = 'point', size=2, aes(group=Category)) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult=1), geom = 'errorbar', size=0.8, aes(group=Category), width=0) +
  theme(axis.text.x = element_text(face= "plain", color= 'black', size=8, angle = 0, hjust = 0.5), axis.text.y = element_text(face = 'plain', color= 'black', size=9)) +
  theme(axis.title.y = element_text(size = 9,face='plain',color='black')) +
  labs(y = "Clathrin spindle recruitment", x = NULL) + 
  ylim(0, 3) +
  theme(legend.position = 'none') +
  theme(legend.title = element_blank()) +
  labs(y = "Clathrin spindle recruitment", x = NULL) +
  expand_limits(y = 0)
plot(combined_clathrin_plot)

combined_gtse1_plot <- ggplot(filter (combined_df), aes(x = Category,y = POI_spindle_ratio, colour = "#ED1C24")) +
  geom_hline(yintercept =1, linetype='dashed', colour='black') +
  geom_quasirandom(alpha = 0.5, stroke = 0) +
  stat_summary(fun.data = mean_se, geom = 'point', size=2, aes(group=Category)) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult=1), geom = 'errorbar', size=0.8, aes(group=Category), width=0) +
  theme(axis.text.x = element_text(face= "plain", color= 'black', size=8, angle = 0, hjust = 0.5), axis.text.y = element_text(face = 'plain', color= 'black', size=9)) +
  theme(axis.title.y = element_text(size = 9,face='plain',color='black')) +
  labs(y = "GTSE1 spindle recruitment", x = NULL) + 
  ylim(0, 3) +
  theme(legend.position = 'none') +
  theme(legend.title = element_blank()) +
  labs(y = "GTSE1 spindle recruitment", x = NULL) +
  expand_limits(y = 0)
plot(combined_gtse1_plot)

# Use cowplot to combine the two plots side by side
combined_plot <- plot_grid(combined_clathrin_plot, combined_gtse1_plot, rel_widths = c(1, 1), rel_heights = c(1,1), ncol = 1)
combined_plot

## Statistics

clathrin_ANOVA <- aov(clathrin_spindle_ratio ~ Category, combined_df)
summary(clathrin_ANOVA)
plot(clathrin_ANOVA,1)
plot(clathrin_ANOVA,2)
TukeyHSD(clathrin_ANOVA, which = "Category")

GTSE1_ANOVA <- aov(POI_spindle_ratio ~ Category, combined_df)
summary(GTSE1_ANOVA)
plot(GTSE1_ANOVA,1)
plot(GTSE1_ANOVA,2)
TukeyHSD(GTSE1_ANOVA, which = "Category")

# Save the plots
ggsave("Output/Plots/combined_plot.pdf", plot = combined_plot, width = 83, height = 132.885, units = 'mm', useDingbats = FALSE)

