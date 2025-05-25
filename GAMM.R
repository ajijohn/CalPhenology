library(dplyr)
library(mgcv)
library(ggplot2)
library(patchwork)
library(svglite)

CrosstabMax <- read.csv('./data/PredictedFlowerAbundanceCrosstab.csv')

reconmeta <- read.csv('./data/reconycx_trailcam_metadata.csv')
reconmetadata  <- reconmeta %>% mutate(camera= paste0("cam",camera) ) %>% as.data.frame()

#Join the tables
FlowerGAM <- CrosstabMax %>% left_join(reconmetadata,by = c("camera"))

#Making a categorical variable out of elevation
FlowerGAM$Elevation <- ifelse(FlowerGAM$elev_m.x > 1600, "high", "low")

# Change "alpine" habitat to "open" and change site names
FlowerGAM$habitat = as.factor(FlowerGAM$habitat.x)
levels(FlowerGAM$habitat)[levels(FlowerGAM$habitat) == "alpine"] <- "open"
levels(FlowerGAM$habitat)

FlowerGAM$site = as.factor(FlowerGAM$site)
levels(FlowerGAM$site)[levels(FlowerGAM$site) == "nbo"] <- "below treeline"
levels(FlowerGAM$site)[levels(FlowerGAM$site) == "subcal"] <- "at treeline"
levels(FlowerGAM$site)[levels(FlowerGAM$site) == "cal"] <- "above treeline"
levels(FlowerGAM$site)

#Remove "no_detections" column
FlowerGAM <- FlowerGAM[, !(names(FlowerGAM) %in% c('no_detections'))]

#Creating a column for elevation habitat combinations
FlowerGAM$Elevation_habitat <- paste(FlowerGAM$site, FlowerGAM$habitat, sep = "/")

# Calculating the daily flower number/camera
FlowerGAM$Flowers <- rowSums(FlowerGAM[, 4:36])

# Calculating the daily species number/camera
FlowerGAM$Species <- apply(FlowerGAM[, 4:36] != 0, 1, sum)


# Fitting the Flower number GAMMs ####
# Fit separate GAMMs for each factor combination
GAMMFlowermodels <- list()

# Loop through each combination and fit a GAMM
unique_combinations <- unique(FlowerGAM$Elevation_habitat)

for (name in unique_combinations) {
  subset_data <- subset(FlowerGAM, Elevation_habitat == name)

  # Fit GAMM for the specific combination with camera as random effect
  model <- gamm(Flowers ~ s(DOY
                            #, bs = "ps", m = c(2,0)# second order smooth, no penalty at zero
                            ),
                random = list(camera = ~1),
                data = subset_data
                ,family = poisson (link = "log")
                )

  # Store the model in the list
  GAMMFlowermodels[[name]] <- model
}

# Summary of models
for (name in unique_combinations) {
  print(paste("Summary for", name))
  print(summary(GAMMFlowermodels[[name]]$gam))  # summary of the smooth terms
}

# Create a sequence of DOY values for prediction
new_data <- data.frame(DOY = seq(min(FlowerGAM$DOY), max(FlowerGAM$DOY), length = 100))

# Initialize an empty data frame to store predictions
predictions <- data.frame()

# Loop through each model and make predictions
for (name in unique_combinations) {
  model <- GAMMFlowermodels[[name]]$gam

  # Predict on link (log) scale
  pred <- predict(model, newdata = new_data, se.fit = TRUE, type = "link")

  # Transform predictions back to response (count) scale
  fit_response <- exp(pred$fit)
  upper <- exp(pred$fit + 1.96 * pred$se.fit)
  lower <- exp(pred$fit - 1.96 * pred$se.fit)

  predictions <- rbind(predictions, data.frame(
    DOY = new_data$DOY,
    Elevation_habitat = name,
    Flowers = fit_response,
    LowerCI = lower,
    UpperCI = upper
  ))
}

#### Fig3c GAMM####
# Define the desired order of levels
desired_order_gam <- c("above treeline, open", "at treeline, canopy", "at treeline, open", "below treeline, canopy", "below treeline, open")

# Map the levels to custom labels
custom_labels_gam <- c("above treeline", "at treeline, canopy", "at treeline, open", "below treeline, canopy", "below treeline, open")

# Define custom labels and colors for the legend
custom_colors <- c( "greenyellow", "gray70", "green2", "gray20", "forestgreen")

# Create the plot
Fig3c <- ggplot() +
  geom_ribbon(data = predictions, aes(x = DOY, ymin = LowerCI, ymax = UpperCI, fill = Elevation_habitat), alpha = 0.3) +
  geom_line(data = predictions, aes(x = DOY, y = Flowers, color = Elevation_habitat)) +
  labs(title = "Flower Abundance",
       x = "Day of Year",
       y = "Flower Number",
       color = "Elevation & Habitat",
       fill = "Elevation & Habitat") +
  scale_color_manual(values = custom_colors, labels = custom_labels_gam) +
  scale_fill_manual(values = custom_colors, labels = custom_labels_gam) +
  ylim(0, NA) +
  theme_minimal() +
  guides(fill = "none", color = guide_legend(title = "Elevation & Habitat"))

Fig3c

ggsave("./figs/Fig3c.svg",Fig3c, width = 25, height = 20, units = "cm")
ggsave("./figs/Fig3c.png",Fig3c, width = 25, height = 20, units = "cm")


# Fitting the Species number GAMMs####
# Fit separate GAMs for each factor combination
GAMMSpeciesModels <- list()

# Loop through each combination and fit a GAM
unique_combinations <- unique(FlowerGAM$Elevation_habitat)

for (name in unique_combinations) {
  subset_data <- subset(FlowerGAM, Elevation_habitat == name)

  model <- gamm(Species ~ s(DOY
                            #, bs = "ps", m = c(2,0)# second order smooth, no penalty at zero
                            ),
                random = list(camera = ~1),
                data = subset_data
                ,family = poisson (link = "log")
                )

  # Store the model in the list
  GAMMSpeciesModels[[name]] <- model
}

# Summary of models
for (name in unique_combinations) {
  print(paste("Summary for", name))
  print(summary(GAMMSpeciesModels[[name]]$gam))  # summary of the smooth terms
}

# Create a sequence of DOY values for prediction
new_data <- data.frame(DOY = seq(min(FlowerGAM$DOY), max(FlowerGAM$DOY), length = 100))

# Initialize an empty data frame to store predictions
predictions <- data.frame()

# Loop through each model and make predictions
for (name in unique_combinations) {
  model <- GAMMSpeciesModels[[name]]$gam

  # Predict on link (log) scale
  pred <- predict(model, newdata = new_data, se.fit = TRUE, type = "link")

  # Transform predictions back to response (count) scale
  fit_response <- exp(pred$fit)
  upper <- exp(pred$fit + 1.96 * pred$se.fit)
  lower <- exp(pred$fit - 1.96 * pred$se.fit)

  predictions <- rbind(predictions, data.frame(
    DOY = new_data$DOY,
    Elevation_habitat = name,
    Species = fit_response,
    LowerCI = lower,
    UpperCI = upper
  ))
}

#Fig3d GAMM####
# Define the desired order of levels
desired_order_gam <- c("above treeline, open", "at treeline, canopy", "at treeline, open", "below treeline, canopy", "below treeline, open")

# Map the levels to custom labels
custom_labels_gam <- c("above treeline", "at treeline, canopy", "at treeline, open", "below treeline, canopy", "below treeline, open")

# Define custom labels and colors for the legend
custom_colors <- c( "greenyellow", "gray70", "green2", "gray20", "forestgreen")

# Create the plot
Fig3d <- ggplot() +
  geom_ribbon(data = predictions, aes(x = DOY, ymin = LowerCI, ymax = UpperCI, fill = Elevation_habitat), alpha = 0.3) +
  geom_line(data = predictions, aes(x = DOY, y = Species, color=Elevation_habitat)) +
  labs(title = "Number of Flowering Species",
       x = "Day of Year",
       y = "Species Number",
       color = "Elevation & Habitat",
       fill = "Elevation & Habitat") +  # Match the legend title for color and fill
  scale_color_manual(values = custom_colors, labels = custom_labels_gam) +
  scale_fill_manual(values = custom_colors, labels = custom_labels_gam) +
  ylim(0, NA) +
  theme_minimal() +
  guides(fill = "none", color = guide_legend(title = "Elevation & Habitat"))  # Hide fill legend, customize color legend


Fig3d
ggsave("./figs/Fig3d.svg",Fig3d, width = 25, height = 20, units = "cm")
ggsave("./figs/Fig3d.png",Fig3d, width = 25, height = 20, units = "cm")

#combine the plots for Fig3 (a & b in ANOVA script)

# Combine the plots into a single layout
combined_plot3 <- (Fig3a | Fig3b) / (Fig3c | Fig3d) + labs(title = "")  +plot_layout(guides = "collect") & plot_annotation(tag_levels = "a")  &
  theme(legend.position = "right")

# Display the combined plot
print(combined_plot3)
ggsave("./figs/Fig3abcd.svg", width = 25, height = 20, units = "cm")
ggsave("./figs/Fig3abcd.png", width = 25, height = 20, units = "cm")
