## GTSE1 fragment data visualisation
## Load in dataframes produced from quantifying individual experiments, then plot graph with stats.
## Original script by James Shelford
## Add-ons and generalisation by Stephen Royle

require(ggplot2)
library(ggbeeswarm)
library(dplyr)
library(multcomp)
library(scales)
library(cowplot)

# make directory for plots if it doesn't exist (it should)
ifelse(!dir.exists("Output"), dir.create("Output"), "Folder exists already")
ifelse(!dir.exists("Output/Plots"), dir.create("Output/Plots"), "Folder exists already")

# Load the dataframes
# these are rds files found in Output/Data
datadir <- "Output/Data"
# for interphase
my_interphase_files <- list.files(datadir,pattern='*interphase.rds',full.names = TRUE)
combined_df_interphase <- bind_rows(lapply(my_interphase_files, readRDS))
combined_df_interphase$Category <- as.factor(combined_df_interphase$Category)

# for mitosis
my_mitosis_files <- list.files(datadir,pattern='*mitosis.rds',full.names = TRUE)
combined_df_mitosis <- bind_rows(lapply(my_mitosis_files, readRDS))
combined_df_mitosis$Category <- as.factor(combined_df_mitosis$Category)

# How many cells in each condition?
summary(combined_df_interphase$Category)
summary(combined_df_mitosis$Category)

# To plot in the correct order, we first need the look-up table
look_up_table <- read.table("Data/lookup.csv", header = TRUE, stringsAsFactors = F, sep = ",")
combined_df_interphase$Category <- factor(combined_df_interphase$Category, levels = look_up_table$Search_category )
combined_df_mitosis$Category <- factor(combined_df_mitosis$Category, levels = look_up_table$Search_category )

# Subsetting the data to remove FKBP-GFP and MAP4
theNames_interphase <- rev(names(combined_df_interphase))
theNames_mitosis <- rev(names(combined_df_mitosis))
theNames_interphase <- theNames_interphase[1:4]
theNames_mitosis <- theNames_mitosis[1:3]

combined_df_interphase <- subset(combined_df_interphase, Category != "FKBP-GFP", select = theNames_interphase)
combined_df_interphase <- subset(combined_df_interphase, Category != "MAP4", select = theNames_interphase)
combined_df_mitosis <- subset(combined_df_mitosis, Category != "FKBP-GFP", select = theNames_mitosis)
combined_df_mitosis <- subset(combined_df_mitosis, Category != "MAP4", select = theNames_mitosis)

## Generating the plots

mitosis_plot <- ggplot(filter(combined_df_mitosis), aes(x=Category, y=POI_spindle_ratio, colour = "#F8766D")) +
  geom_hline(yintercept =1, linetype='dashed', colour='black') +
  geom_quasirandom(alpha=0.5, stroke=0) + 
  stat_summary(fun.data = mean_se, geom = 'point', size=2, aes(group=Category))+
  stat_summary(fun.data = mean_sdl, fun.args = list(mult=1), geom = 'errorbar', size=0.8, aes(group=Category), width=0) +
  theme(axis.text.x = element_blank(),axis.text.y = element_text(face = 'plain', color= 'black', size=9)) +
  theme(axis.title.y = element_text(size = 9,face='plain',color='black')) +
  scale_colour_manual(values = "#F8766D") + 
  labs(y = "GTSE1 spindle recruitment", x = NULL) + 
  theme(legend.position = 'none') +
  theme(legend.title = element_blank()) +
  ylim(0,3)

interphase_plot <- ggplot(data = combined_df_interphase, aes(x = Category,y = POI_MT_ratio, colour = '#00BFC4')) +
  geom_quasirandom(alpha = 0.5, stroke = 0) +
  geom_hline(yintercept =1, linetype='dashed', colour='black') +
  stat_summary(fun.data = mean_se, geom = 'point', size=2, aes(group=Category))+
  stat_summary(fun.data = mean_sdl, fun.args = list(mult=1), geom = 'errorbar', size=0.8, aes(group=Category), width=0) +
  theme(axis.text.x = element_text(face= "plain", color= 'black', size=8.5, angle = 40, hjust = 0.5, vjust = 0.6), axis.text.y = element_text(face = 'plain', color= 'black', size=9)) +
  theme(axis.title.y = element_text(size = 9,face='plain',color='black')) +
  scale_color_manual(values = '#00BFC4') +
  labs(y = "Ratio (intensity MT / cytoplasm)", x = NULL) +
  theme(legend.position = 'none') +
  theme(legend.title = element_blank()) +
  ylim(0,7)

# combine plots
combined_plot <- plot_grid(mitosis_plot, interphase_plot, ncol = 1)

## Statistics
interphase_ANOVA <- aov(POI_MT_ratio ~ Category, combined_df_interphase)
summary(interphase_ANOVA)
plot(interphase_ANOVA,1)
plot(interphase_ANOVA,2)
TukeyHSD(interphase_ANOVA, which = "Category")

mitosis_ANOVA <- aov(POI_spindle_ratio ~ Category, combined_df_mitosis)
summary(mitosis_ANOVA)
plot(mitosis_ANOVA,1)
plot(mitosis_ANOVA,2)
TukeyHSD(mitosis_ANOVA, which = "Category")

# Save the plots
ggsave("Output/Plots/mitosis_plot.png", plot = mitosis_plot, dpi = 300)
ggsave("Output/Plots/interphase_plot.png", plot = interphase_plot, dpi = 300)
ggsave("Output/Plots/GTSE1_fragment_combined.pdf", plot = combined_plot, width = 78, height = 140, units = 'mm', useDingbats = FALSE)

