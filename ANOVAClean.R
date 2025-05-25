#This script performs ANOVAs and plots the results for the CalPhe 2022 flowering phenology data
library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(ggsignif)
library(patchwork)
library(svglite)
library(emmeans)
library(multcomp)
#Get the data
EarlyLateSeason <- read.csv('./data/EarlyLateSeason.csv')

#Alternative data set without extrapolated flowering in the end and beginning of the observation period
#EarlyLateSeason <- read.csv('C:/Users/Mikko/Desktop/CalPhe/2022/CalPheR/FinalR/Data/EarlyLateSeasonNoTails.csv')

# Change "alpine" habitat to "open" and sites to treeline-relevant labels
EarlyLateSeason$habitat = as.factor(EarlyLateSeason$habitat)
levels(EarlyLateSeason$habitat)[levels(EarlyLateSeason$habitat) == "alpine"] <- "open"
levels(EarlyLateSeason$habitat)
EarlyLateSeason$site[EarlyLateSeason$site == "nbo"] <- "below treeline"
EarlyLateSeason$site[EarlyLateSeason$site == "subcal"] <- "at treeline"
EarlyLateSeason$site[EarlyLateSeason$site == "cal"] <- "above treeline"

EarlyLateSeason$SiteHabitat <- paste(EarlyLateSeason$site, EarlyLateSeason$habitat, sep = "/")


#ANOVAs####
# Early ANOVA
EarlyLateSeason %>%
  group_by(SiteHabitat) %>%
  summarise(
    AverageEarly = format(mean(Early, na.rm = TRUE), digits = 5, nsmall = 5),
    StdDevEarly = format(sd(Early, na.rm = TRUE), digits = 5, nsmall = 5)
  )

ANOVAEarly <- aov(Early ~ SiteHabitat, data = EarlyLateSeason)
summary(ANOVAEarly)
# non-significant

# Late ANOVA
EarlyLateSeason %>%
  group_by(SiteHabitat) %>%
  summarise(
    AverageEarly = format(mean(Late, na.rm = TRUE), digits = 5, nsmall = 5),
    StdDevEarly = format(sd(Late, na.rm = TRUE), digits = 5, nsmall = 5)
  )

ANOVALate <- aov(Late ~ SiteHabitat, data = EarlyLateSeason)
summary(ANOVALate)
# significant

# Average flowering ANOVA
EarlyLateSeason %>%
  group_by(SiteHabitat) %>%
  summarise(
    AverageEarly = format(mean(AverageFlowering, na.rm = TRUE), digits = 5, nsmall = 5),
    StdDevEarly = format(sd(AverageFlowering, na.rm = TRUE), digits = 5, nsmall = 5)
  )

ANOVAAverageFlowering <- aov(AverageFlowering ~ SiteHabitat, data = EarlyLateSeason)
summary(ANOVAAverageFlowering)
# non-significant

# Season length ANOVA
EarlyLateSeason %>%
  group_by(SiteHabitat) %>%
  summarise(
    AverageEarly = format(mean(SeasonLength, na.rm = TRUE), digits = 5, nsmall = 5),
    StdDevEarly = format(sd(SeasonLength, na.rm = TRUE), digits = 5, nsmall = 5)
  )

ANOVASeasonLength <- aov(SeasonLength ~ SiteHabitat, data = EarlyLateSeason)
summary(ANOVASeasonLength)
# significant

# Flower Number ANOVA
EarlyLateSeason %>%
  group_by(SiteHabitat) %>%
  summarise(
    AverageEarly = format(mean(FlowerNumber, na.rm = TRUE), digits = 5, nsmall = 5),
    StdDevEarly = format(sd(FlowerNumber, na.rm = TRUE), digits = 5, nsmall = 5)
  )

ANOVAFlowerNumber <- aov(FlowerNumber ~ SiteHabitat, data = EarlyLateSeason)
summary(ANOVAFlowerNumber)
# significant

# Species Number ANOVA
EarlyLateSeason %>%
  group_by(SiteHabitat) %>%
  summarise(
    AverageEarly = format(mean(SpeciesNumber, na.rm = TRUE), digits = 5, nsmall = 5),
    StdDevEarly = format(sd(SpeciesNumber, na.rm = TRUE), digits = 5, nsmall = 5)
  )

ANOVASpeciesNumber <- aov(SpeciesNumber ~ SiteHabitat, data = EarlyLateSeason)
summary(ANOVASpeciesNumber)
# significant


#Post hoc tests####
#Late
# Perform Tukey's HSD post hoc test
PosthocLate <- TukeyHSD(ANOVALate)
print(PosthocLate)
plot(PosthocLate)

#Season length
# Perform Tukey's HSD post hoc test
PosthocSeasonLength <- TukeyHSD(ANOVASeasonLength)
print(PosthocSeasonLength)
plot(PosthocSeasonLength)

#Flower Number
# Perform Tukey's HSD post hoc test
PosthocFlowerNumber <- TukeyHSD(ANOVAFlowerNumber)
print(PosthocFlowerNumber)
plot(PosthocFlowerNumber)

#Species Number
# Perform Tukey's HSD post hoc test
PosthocSpeciesNumber <- TukeyHSD(ANOVASpeciesNumber)
print(PosthocSpeciesNumber)
plot(PosthocSpeciesNumber)


#Plots####
#Fig 2####
# Early
# Reorder SiteHabitat levels
EarlyLateSeason$SiteHabitat <- factor(
  EarlyLateSeason$SiteHabitat,
  levels = c(
    "below treeline/open",
    "below treeline/canopy",
    "at treeline/open",
    "at treeline/canopy",
    "above treeline/open"
  )
)

# Ensure Habitat is a factor with the correct levels
EarlyLateSeason$Site <- factor(
  EarlyLateSeason$Site,
  levels = c("below treeline", "at treeline", "above treeline")
)

# Create the plot with customized colors
Fig2a <- ggplot(EarlyLateSeason, aes(x = SiteHabitat, y = Early, fill = habitat)) +
  geom_boxplot(outlier.shape = NA) +
  #stat_compare_means(method = "anova") +  # Adds overall ANOVA p-value
  geom_jitter(width = 0.2,
              height = 0,
              alpha = 1,
              color = "black"
  ) +
  scale_fill_manual(values = c("open" = "palegreen2", "canopy" = "gray60")) +  # Custom colors
  theme_minimal() +
  labs(
    title = "Earliest flowering of species",
    x = "Study site",
    y = "Day of the Year"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+  # Rotate x-axis labels for better visibility
  theme(legend.position="none")

Fig2a

# Late
# Perform Tukey's HSD post hoc test
tukey_data <- as.data.frame(PosthocLate$SiteHabitat)

# Extract significant pairwise comparisons
sig_pairs <- subset(tukey_data, `p adj` < 0.05)
sig_pairs$groups <- rownames(sig_pairs)

# Split groups into separate columns for comparisons
sig_pairs <- sig_pairs %>%
  separate(groups, into = c("group1", "group2"), sep = "-")

# Prepare the list of comparisons
comparisons <- split(sig_pairs[, c("group1", "group2")], seq(nrow(sig_pairs)))
comparisons <- lapply(comparisons, as.character)

# Reorder SiteHabitat levels
EarlyLateSeason$SiteHabitat <- factor(
  EarlyLateSeason$SiteHabitat,
  levels = c(
    "below treeline/open",
    "below treeline/canopy",
    "at treeline/open",
    "at treeline/canopy",
    "above treeline/open"
  )
)

# Ensure Habitat is a factor with the correct levels
EarlyLateSeason$Site <- factor(
  EarlyLateSeason$Site,
  levels = c("below treeline", "at treeline", "above treeline")
)

# Create the plot with customized colors
Fig2b <- ggplot(EarlyLateSeason, aes(x = SiteHabitat, y = Late, fill = habitat)) +
  geom_boxplot(outlier.shape = NA) +
  #stat_compare_means(method = "anova") +  # Adds overall ANOVA p-value
  geom_signif(
    comparisons = comparisons,  # Pass the list of comparisons
    map_signif_level = TRUE,
    textsize = 3
  ) +
  geom_jitter(width = 0.2,
              height = 0,
              alpha = 1,
              color = "black"
  ) +
  scale_fill_manual(values = c("open" = "palegreen2", "canopy" = "gray60")) +  # Custom colors
  theme_minimal() +
  labs(
    title = "Latest flowering of species",
    x = "Study site",
    y = "Day of the Year",
    fill = "Habitat"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better visibility

Fig2b

# Average flowering
# Reorder SiteHabitat levels
EarlyLateSeason$SiteHabitat <- factor(
  EarlyLateSeason$SiteHabitat,
  levels = c(
    "below treeline/open",
    "below treeline/canopy",
    "at treeline/open",
    "at treeline/canopy",
    "above treeline/open"
  )
)

# Ensure Habitat is a factor with the correct levels
EarlyLateSeason$Site <- factor(
  EarlyLateSeason$Site,
  levels = c("below treeline", "at treeline", "above treeline")
)

# Create the plot with customized colors
Fig2c <- ggplot(EarlyLateSeason, aes(x = SiteHabitat, y = AverageFlowering, fill = habitat)) +
  geom_boxplot(outlier.shape = NA) +
  #stat_compare_means(method = "anova") +  # Adds overall ANOVA p-value
  geom_jitter(width = 0.2,
              height = 0,
              alpha = 1,
              color = "black"
  ) +
  scale_fill_manual(values = c("open" = "palegreen2", "canopy" = "gray60")) +  # Custom colors
  theme_minimal() +
  labs(
    title = "Average flowering of the community",
    x = "Study site",
    y = "Day of the Year"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+  # Rotate x-axis labels for better visibility
  theme(legend.position="none")

Fig2c

# Season Length
# Perform Tukey's HSD post hoc test
tukey_data <- as.data.frame(PosthocSeasonLength$SiteHabitat)

# Extract significant pairwise comparisons
sig_pairs <- subset(tukey_data, `p adj` < 0.05)
sig_pairs$groups <- rownames(sig_pairs)

# Split groups into separate columns for comparisons
sig_pairs <- sig_pairs %>%
  separate(groups, into = c("group1", "group2"), sep = "-")

# Prepare the list of comparisons
comparisons <- split(sig_pairs[, c("group1", "group2")], seq(nrow(sig_pairs)))
comparisons <- lapply(comparisons, as.character)

# Reorder SiteHabitat levels
EarlyLateSeason$SiteHabitat <- factor(
  EarlyLateSeason$SiteHabitat,
  levels = c(
    "below treeline/open",
    "below treeline/canopy",
    "at treeline/open",
    "at treeline/canopy",
    "above treeline/open"
  )
)

# Ensure Habitat is a factor with the correct levels
EarlyLateSeason$site <- factor(
  EarlyLateSeason$Site,
  levels = c("below treeline", "at treeline", "above treeline")
)

# Create the plot with customized colors
Fig2d <- ggplot(EarlyLateSeason, aes(x = SiteHabitat, y = SeasonLength, fill = habitat)) +
  geom_boxplot(outlier.shape = NA) +
  #stat_compare_means(method = "anova") +  # Adds overall ANOVA p-value
  geom_signif(
    comparisons = comparisons,  # Pass the list of comparisons
    map_signif_level = TRUE,
    textsize = 3,
    y_position = c(130,140) #Adjust y_position to avoid overlap
  ) +
  geom_jitter(width = 0.2,
              height = 0,
              alpha = 1,
              color = "black"
  ) +
  scale_fill_manual(values = c("open" = "palegreen2", "canopy" = "gray60")) +  # Custom colors
  theme_minimal() +
  labs(
    title = "Flowering season length",
    x = "Study site",
    y = "Days"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+ theme(legend.position="none")  # Rotate x-axis labels for better visibility

Fig2d

# Combine the plots using patchwork
combined_plot2 <- (Fig2a | Fig2b) / (Fig2c | Fig2d)+ labs(title = "") +   plot_layout(guides = "collect") & plot_annotation(tag_levels = "a")
combined_plot2
ggsave("./figs/combined_plot2.svg",combined_plot2, width = 25, height = 20, units = "cm")
ggsave("./figs/combined_plot2.png",combined_plot2, width = 25, height = 20, units = "cm")

#Fig 3 a b####

# Flower number
# Perform Tukey's HSD post hoc test
tukey_data <- as.data.frame(PosthocFlowerNumber$SiteHabitat)

# Extract significant pairwise comparisons
sig_pairs <- subset(tukey_data, `p adj` < 0.05)
sig_pairs$groups <- rownames(sig_pairs)

# Split groups into separate columns for comparisons
sig_pairs <- sig_pairs %>%
  separate(groups, into = c("group1", "group2"), sep = "-")

# Prepare the list of comparisons
comparisons <- split(sig_pairs[, c("group1", "group2")], seq(nrow(sig_pairs)))
comparisons <- lapply(comparisons, as.character)

# Reorder SiteHabitat levels
EarlyLateSeason$SiteHabitat <- factor(
  EarlyLateSeason$SiteHabitat,
  levels = c(
    "below treeline/open",
    "below treeline/canopy",
    "at treeline/open",
    "at treeline/canopy",
    "above treeline/open"
  )
)

# Ensure Habitat is a factor with the correct levels
EarlyLateSeason$Site <- factor(
  EarlyLateSeason$Site,
  levels = c("below treeline", "at treeline", "above treeline")
)

# Create the plot with customized colors
Fig3a <- ggplot(EarlyLateSeason, aes(x = SiteHabitat, y = FlowerNumber, fill = habitat)) +
  geom_boxplot(outlier.shape = NA) +
  #stat_compare_means(method = "anova") +  # Adds overall ANOVA p-value
  geom_signif(
    comparisons = comparisons,  # Pass the list of comparisons
    map_signif_level = TRUE,
    textsize = 3,
    y_position = c(1700) #Adjust y_position to avoid overlap
  ) +
  geom_jitter(width = 0.2,
              height = 0,
              alpha = 1,
              color = "black"
  ) +
  scale_fill_manual(values = c("open" = "palegreen2", "canopy" = "gray60")) +  # Custom colors
  theme_minimal() +
  labs(
    title = "Flower Number",
    x = "Study site",
    y = "Flower Number",
    fill = "Habitat"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better visibility
Fig3a

# Species number
# Perform Tukey's HSD post hoc test
tukey_data <- as.data.frame(PosthocSpeciesNumber$SiteHabitat)

# Extract significant pairwise comparisons
sig_pairs <- subset(tukey_data, `p adj` < 0.05)
sig_pairs$groups <- rownames(sig_pairs)

# Split groups into separate columns for comparisons
sig_pairs <- sig_pairs %>%
  separate(groups, into = c("group1", "group2"), sep = "-")

# Prepare the list of comparisons
comparisons <- split(sig_pairs[, c("group1", "group2")], seq(nrow(sig_pairs)))
comparisons <- lapply(comparisons, as.character)

# Reorder SiteHabitat levels
EarlyLateSeason$SiteHabitat <- factor(
  EarlyLateSeason$SiteHabitat,
  levels = c(
    "below treeline/open",
    "below treeline/canopy",
    "at treeline/open",
    "at treeline/canopy",
    "above treeline/open"
  )
)

# Ensure Habitat is a factor with the correct levels
EarlyLateSeason$Site <- factor(
  EarlyLateSeason$Site,
  levels = c("below treeline", "at treeline", "above treeline")
)

# Create the plot with customized colors
Fig3b <- ggplot(EarlyLateSeason, aes(x = SiteHabitat, y = SpeciesNumber, fill = habitat)) +
  geom_boxplot(outlier.shape = NA) +
  #stat_compare_means(method = "anova") +  # Adds overall ANOVA p-value
  geom_signif(
    comparisons = comparisons,  # Pass the list of comparisons
    map_signif_level = TRUE,
    textsize = 3,
    y_position = c(22,21,19,20) #Adjust y_position to avoid overlap
  ) +
  geom_jitter(width = 0.2,
              height = 0,
              alpha = 1,
              color = "black"
  ) +
  scale_fill_manual(values = c("open" = "palegreen2", "canopy" = "gray60")) +  # Custom colors
  theme_minimal() +
  labs(
    title = "Species Number",
    x = "Study site",
    y = "Species Number",
    fill = "Habitat"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better visibility
Fig3b

# Combine the plots using patchwork
combined_plot3ab <- (Fig3a | Fig3b) + labs(title = "")  & plot_annotation(tag_levels = "a")
combined_plot3ab
ggsave("./figs/combined_plot3ab.svg",combined_plot3ab, width = 25, height = 15, units = "cm")
ggsave("./figs/combined_plot3ab.png",combined_plot3ab, width = 25, height = 15, units = "cm")

