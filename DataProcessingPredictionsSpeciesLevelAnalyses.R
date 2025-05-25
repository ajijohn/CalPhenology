# This script compiles, and modifies the CalPhe 2022 Flower count data
# Also analyses and plots the species level patterns in flowering phenology
library(tidyverse)
library(lubridate)
library(tidyr)
library(dplyr)
library(ggplot2)
library(zoo)
library(reshape2)
library(patchwork)
# Set the working directory to the folder containing the prediction CSV files
#setwd("./data")

# Set the directory for the data
data_dir <- "./data/"

# Set the directory containing CSV files
csv_dir <- "./data/Predictions-Colab-last-weights/"

# List all CSV files in the directory (full paths)
csv_files <- list.files(path = csv_dir, pattern = "\\.csv$", full.names = TRUE)
# Initialize an empty dataframe
AllPredictions <- data.frame()

# Read and combine all CSV files
for (file in csv_files) {
  temp_df <- read.csv(file, stringsAsFactors = FALSE)
  AllPredictions <- rbind(AllPredictions, temp_df)
}
# Creating a Cross tab with species as columns ----------------------------

#All image metadata
ct <- read_csv('./data/CTAllNonfilteredmeta.csv')

#Flower predictions
predictions <- AllPredictions

predictionsbeefed <- predictions %>%
  separate(file, c("base","imgfolder"),
           "/content/drive/MyDrive/Calanda-Drive/Trail Cam Photos"
           )%>%

  separate(imgfolder,c("folder","subfolder","file","filef"),"/")

names(predictionsbeefed) <- c( "mbase",
                               "mfolder",
                               "camera",
                               "subfolder",
                               "filename",
                               "img",
                               "class",
                               "confidence")

head(predictionsbeefed)

# Camera metadata
reconmeta <- read_csv('./data/reconycx_trailcam_metadata.csv')
reconmetadata  <- reconmeta %>% mutate(camera= paste0("cam",camera) ) %>% as.data.frame()

predandcamerameta <- predictionsbeefed %>% left_join(reconmetadata,by = c("camera"))

predandcamerametafilemeta <-ct %>% mutate(filename=paste0(filename,".JPG")) %>%
  inner_join(predandcamerameta)
head(predandcamerametafilemeta)

# Creating the crosstab
Crosstab <- predandcamerametafilemeta %>%
  group_by(camera,
           filedate,
           filetime,
           elev_m,
           class) %>%
  summarise(count = n())%>%
  spread(key = class, value = count, fill = 0)

# Change the dates into Day of the Year
# Convert "DD-MM-YYYY" to Date class
Crosstab$filedate <- dmy(Crosstab$filedate)

# Convert Date to Day of the Year (DOY)
Crosstab$DOY <- as.numeric(Crosstab$filedate - floor_date(Crosstab$filedate, "year")) + 1

#Replace spaces in column names with underscores
colnames(Crosstab) <- gsub(" ", "_", colnames(Crosstab))

write.csv(Crosstab,'./data/CrosstabRAW.csv')

# Remove lonely flower observations (likely false positive annotations)
# This process with a filter of +-5 days removes 5 species from the data
Crosstab <- read.csv('./data/CrosstabRAW.csv')
#Remove "no_detections" column
Crosstab <- Crosstab[, !(names(Crosstab) %in% c('X'))]

# Function to check if there are at least two non-zero observations in a column
check_non_zero <- function(column_values) {
  non_zero_count <- sum(column_values != 0, na.rm = TRUE)
  return(non_zero_count >= 2)
}

# Get unique values of camera column
unique_cameras <- unique(Crosstab$camera)

# Iterate over unique camera values
for (camera in unique_cameras) {
  # Subset data for the current camera
  camera_data <- Crosstab[Crosstab$camera == camera, ]

  # Iterate over columns from 6 to 43
  for (col in 5:42) {
    # Iterate over rows
    for (i in 1:nrow(camera_data)) {
      # Get DOY value for the current row
      current_doy <- camera_data[i, "DOY"]

      # Extract subset of rows within the range of -5 to 5 from current DOY
      subset_rows <- camera_data[camera_data$DOY >= (current_doy - 5) &
                                   camera_data$DOY <= (current_doy + 5), col]

      # Check if there are at least two non-zero observations in the subset
      if (!check_non_zero(subset_rows)) {
        camera_data[i, col] <- 0
      }
    }
  }

  # Update the original data frame with modified data for the current camera
  Crosstab[Crosstab$camera == camera, 5:42] <- camera_data[, 5:42]
}

write.csv(Crosstab,'./data/Crosstab0.csv')

# Only keeping the daily max for all the species
CrosstabMax <- Crosstab %>%
  group_by(camera, filedate, elev_m) %>%
  summarize(across(starts_with("Achillea_millefolium"):last_col(), ~ max(.)))

#Remove "no_detections" column
CrosstabMax <- CrosstabMax[, !(names(CrosstabMax) %in% c('no_detections'))]

# Remove species columns without data
# Check which columns have all zero values
all_zero_cols <- sapply(CrosstabMax, function(x) all(x == 0))
# Identify columns with at least one non-zero value
non_zero_cols <- !all_zero_cols
# Subset CrosstabMax to keep only columns with at least one non-zero value
CrosstabMax <- CrosstabMax[, non_zero_cols]

write.csv(CrosstabMax,'./data/CrosstabMax.csv')


# Estimating and predicting flowering curves ---------------------------------------------
# Here we interpolate and extrapolate local flower numbers based on the observations

# First, extrapolation of the values that are either earlier or later than the observation period
# Earliest DOY for images at each camera
CrosstabRaw <- read.csv('./data/CrosstabRaw.csv')

# Find the earliest DOY value for each unique value in column "camera"
EarliestPhotos <- aggregate(DOY ~ camera, data = CrosstabRaw, FUN = min)

# Find the latest DOY value for each unique value in column "camera"
LatestPhotos <- aggregate(DOY ~ camera, data = CrosstabRaw, FUN = max)

# The earliest and latest flowering for each species
CrosstabMax <- read.csv('./data/CrosstabMax.csv')
#Remove "X" column
CrosstabMax <- CrosstabMax[, !(names(CrosstabMax) %in% c('X'))]


# Identifying the earliest DOY for all the species
species_columns <- colnames(CrosstabMax[4:36])
# Create an empty data frame to store the results
CrosstabEarly <- unique(CrosstabMax[, c("camera", "elev_m")])

# Iterate over each species
for (species in species_columns) {
  # Find the smallest DOY for the current species for each camera and elev_m combination
  min_doy_species <- aggregate(DOY ~ camera + elev_m, data = subset(CrosstabMax, CrosstabMax[[species]] > 0), FUN = min)

  # Merge the result with the original data frame using a left join
  CrosstabEarly <- merge(CrosstabEarly, min_doy_species, by.x = c("camera", "elev_m"), by.y = c("camera", "elev_m"), all.x = TRUE)

  # Rename the new column with the original species name
  colnames(CrosstabEarly)[ncol(CrosstabEarly)] <- species
}


# Identifying the latest DOY for all the species
# Create an empty data frame to store the results
CrosstabLate <- unique(CrosstabMax[, c("camera", "elev_m")])
# Iterate over each species
for (species in species_columns) {
  # Find the smallest DOY for the current species for each camera and elev_m combination
  max_doy_species <- aggregate(DOY ~ camera + elev_m, data = subset(CrosstabMax, CrosstabMax[[species]] > 0), FUN = max)

  # Merge the result with the original data frame using a left join
  CrosstabLate <- merge(CrosstabLate, max_doy_species, by.x = c("camera", "elev_m"), by.y = c("camera", "elev_m"), all.x = TRUE)

  # Rename the new column with the original species name
  colnames(CrosstabLate)[ncol(CrosstabLate)] <- species
}


head(CrosstabEarly)
head(CrosstabLate)


# Calculate the difference between species earliest observation and the earliest photo
relevant_columns <- CrosstabEarly[, 3:35]
# Loop through each column and subtract the values from EarliestPhotos$DOY
for (i in 1:ncol(relevant_columns)) {
  relevant_columns[, i] <- relevant_columns[, i] - EarliestPhotos$DOY
}

FloweringStart <- cbind(EarliestPhotos, relevant_columns)

# Pivot FloweringStart data frame
FloweringStartPivot <- pivot_longer(
  data = FloweringStart,
  cols = 3:35,
  names_to = "Species",
  values_to = "FloweringPhotoStart"
)


# Calculate the difference between species latest observation and the latest photo
relevant_columns <- CrosstabLate[, 3:35]
# Loop through each column and subtract the values from EarliestPhotos$DOY
for (i in 1:ncol(relevant_columns)) {
  relevant_columns[, i] <- LatestPhotos$DOY - relevant_columns[, i]
}

FloweringEnd <- cbind(LatestPhotos, relevant_columns)

# Pivot FloweringEnd data frame
FloweringEndPivot <- pivot_longer(
  data = FloweringEnd,
  cols = 3:35,
  names_to = "Species",
  values_to = "FloweringPhotoEnd"
)

# Find out the slope of flower abundance as a function of DOY for the earliest flower observations
# Do this for species that were present in the very early season

CrosstabMax <- read.csv('./data/CrosstabMax.csv')
#Remove "X" column
CrosstabMax <- CrosstabMax[, !(names(CrosstabMax) %in% c('X'))]

# Get unique values in column "camera"
unique_cameras <- unique(CrosstabMax$camera)

# Initialize an empty dataframe to store the results
result_df <- data.frame(matrix(ncol = 12, nrow = 0))
colnames(result_df) <- c("camera", "Species", "value1", "DOY1", "value2", "DOY2", "value3", "DOY3", "value4", "DOY4", "AbundanceEarly", "AbundanceLate")

# Loop through each unique camera
for (camera in unique_cameras) {
  # Filter the dataframe to include only rows where camera is the current camera
  camera_data <- CrosstabMax[CrosstabMax$camera == camera, ]

  # Loop through columns 4 to 36
  for (col in 4:36) {
    # Filter the dataframe to include only rows where the current column is non-zero
    non_zero_data <- camera_data[camera_data[, col] != 0, ]

    # Sort the filtered dataframe by the column "DOY"
    sorted_data <- non_zero_data[order(non_zero_data$DOY), ]

    # Extract the two smallest non-zero values and their respective DOY values
    smallest_values <- sorted_data[1:2, c(names(sorted_data)[col], "DOY")]

    # Extract the two largest non-zero values and their respective DOY values
    largest_values <- tail(sorted_data[sorted_data[, col] != 0, c(names(sorted_data)[col], "DOY")], 2)

    # Calculate the sum of non-zero values
    sum_values <- sum(sorted_data[, col])

    # Extract the column name of the original dataframe
    column_name <- names(sorted_data)[col]

    # Create a dataframe with the results
    result <- data.frame(
      camera = camera,
      Species = column_name,
      value1 = smallest_values[1, 1],
      DOY1 = smallest_values[1, 2],
      value2 = ifelse(nrow(smallest_values) > 1, smallest_values[2, 1], NA),
      DOY2 = ifelse(nrow(smallest_values) > 1, smallest_values[2, 2], NA),
      value3 = largest_values[1, 1],
      DOY3 = largest_values[1, 2],
      value4 = largest_values[2, 1],
      DOY4 = largest_values[2, 2],
      AbundanceEarly = smallest_values[1, 1],
      AbundanceLate = largest_values[2, 1]
    )

    # Bind the result to the main dataframe
    result_df <- rbind(result_df, result)
  }
}

# Sort the result dataframe by "camera" and "Species"
SpeciesSlope <- result_df[order(result_df$camera, result_df$Species), ]

# Calculate the slope for each row
SpeciesSlope$EarlySlope <- (SpeciesSlope$value2 - SpeciesSlope$value1) / (SpeciesSlope$DOY2 - SpeciesSlope$DOY1)
SpeciesSlope$LateSlope <- (SpeciesSlope$value4 - SpeciesSlope$value3) / (SpeciesSlope$DOY4 - SpeciesSlope$DOY3)

# Standardise the slopes by dividing with the local species abundance
SpeciesSlope$EarlySlopeSTD <- SpeciesSlope$EarlySlope / SpeciesSlope$AbundanceEarly
SpeciesSlope$LateSlopeSTD <- SpeciesSlope$LateSlope / SpeciesSlope$AbundanceLate

# Estimate a Gaussian distribution along which the extrapolated values would be located
# Define parameters
mu <- 5     # Mean of the Gaussian curve
sigma <- 2  # Standard deviation of the Gaussian curve
x <- seq(0, 10, by = 0.1)  # Range of x values

# Calculate the Gaussian curve values
gaussian_curve <- dnorm(x, mean = mu, sd = sigma)

# Calculate the derivative of the Gaussian curve
dx <- diff(x[1:2])
gaussian_derivative <- diff(gaussian_curve) / dx

# Add a value of 0 to position 51
gaussian_derivative <- c(gaussian_derivative[1:50], 0, gaussian_derivative[51:100])

# Extract values for integer x-values
integer_x_values <- x[x %% 1 == 0]
integer_gaussian_curve <- gaussian_curve[x %% 1 == 0]
integer_gaussian_derivative <- gaussian_derivative[x %% 1 == 0]

# Create a data frame to store the density and derivative values for integer x-values
GaussianCurve <- data.frame(x = integer_x_values,
                            density = integer_gaussian_curve,
                            derivative = integer_gaussian_derivative)

# Plot both the Gaussian curve and its derivative
par(mfrow = c(2, 1), mar = c(4, 4, 2, 2))  # Set up multiple plots
plot(x, gaussian_curve, type = "l", col = "blue", xlab = "x", ylab = "Density", main = "Gaussian Curve")
points(integer_x_values, integer_gaussian_curve, col = "red", pch = 16)  # Add points for integer x-values

plot(x[-1], gaussian_derivative[-1], type = "l", col = "blue", xlab = "x", ylab = "Derivative", main = "Derivative of Gaussian Curve")
points(integer_x_values[-1], integer_gaussian_derivative[-1], col = "red", pch = 16)  # Add points for integer x-values

# Remove tail values of the derivative curve (to avoid multiple peaks)
positions_to_replace <- c(1, 2, 3, 9, 10, 11)
GaussianCurve$derivative[positions_to_replace] <- NA


# Find the position of the flowering curve slope on a Gaussian curve
# Initialize GaussianDOY column in SpeciesSlope dataframe
SpeciesSlope$GaussianDOYEarly <- NA

# Iterate over each row in SpeciesSlope
for (i in 1:nrow(SpeciesSlope)) {
  # Get the value of EarlySlopeSTD in the current row
  slope_std <- SpeciesSlope$EarlySlopeSTD[i]

  # Find the index of the closest numeric value in GaussianCurve$derivative
  closest_index <- which.min(abs(GaussianCurve$derivative - slope_std))

  # Check if closest_index is not empty
  if (length(closest_index) > 0) {
    # Get the corresponding x value from GaussianCurve$x
    closest_x <- GaussianCurve$x[closest_index]

    # Write the closest x value into the GaussianDOY column in SpeciesSlope
    SpeciesSlope$GaussianDOYEarly[i] <- closest_x
  } else {
    # If no matching index is found, set GaussianDOY to NA
    SpeciesSlope$GaussianDOYEarly[i] <- NA
  }
}


# Find the position of the flowering curve slope on a Gaussian curve
# Initialize GaussianDOY column in SpeciesSlope dataframe
SpeciesSlope$GaussianDOYLate <- NA

# Iterate over each row in SpeciesSlope
for (i in 1:nrow(SpeciesSlope)) {
  # Get the value of EarlySlopeSTD in the current row
  slope_std <- SpeciesSlope$LateSlopeSTD[i]

  # Find the index of the closest numeric value in GaussianCurve$derivative
  closest_index <- which.min(abs(GaussianCurve$derivative - slope_std))

  # Check if closest_index is not empty
  if (length(closest_index) > 0) {
    # Get the corresponding x value from GaussianCurve$x
    closest_x <- GaussianCurve$x[closest_index]

    # Write the closest x value into the GaussianDOY column in SpeciesSlope
    SpeciesSlope$GaussianDOYLate[i] <- closest_x
  } else {
    # If no matching index is found, set GaussianDOY to NA
    SpeciesSlope$GaussianDOYLate[i] <- NA
  }
}

# Fill in the number of values for earlier DOYs
# Create a sequence of integers from 0 to 10
sequence <- 1:10

# Generate column names
new_column_names <- paste0("DOYE", sequence)

# Create empty columns in SpeciesSlope data frame
SpeciesSlope[, new_column_names] <- NA

# Loop through each row and fill in the new columns
for (i in 1:nrow(SpeciesSlope)) {
  # Get the value of GaussianDOYEarly for the current row
  gaussian_doy_early <- SpeciesSlope$GaussianDOYEarly[i]

  # Check if GaussianDOYEarly is a finite number
  if (is.finite(gaussian_doy_early)) {
    # Create a sequence of integers smaller than GaussianDOYEarly, but reverse the sequence
    smaller_values <- rev(seq(0, gaussian_doy_early - 1, by = 1))

    # Assign the smaller values to the corresponding columns
    SpeciesSlope[i, new_column_names[1:length(smaller_values)]] <- smaller_values
  } else {
    # If GaussianDOYEarly is not a finite number, fill the corresponding columns with NA
    SpeciesSlope[i, new_column_names] <- NA
  }
}

# Fill in the number of values for later DOYs
# Create a sequence of integers from 1 to 10
sequence <- 1:10

# Generate column names
new_column_names <- paste0("DOYL", sequence)

# Create empty columns in SpeciesSlope data frame
SpeciesSlope[, new_column_names] <- NA

# Loop through each row and fill in the new columns
for (i in 1:nrow(SpeciesSlope)) {
  # Get the value of GaussianDOYLate for the current row
  gaussian_doy_late <- SpeciesSlope$GaussianDOYLate[i]

  # Check if GaussianDOYLate is a finite number
  if (is.finite(gaussian_doy_late)) {
    # Create a sequence of integers larger than GaussianDOYLate, up to 10
    larger_values <- seq(gaussian_doy_late + 1, 10, by = 1)

    # Assign the larger values to the corresponding columns
    SpeciesSlope[i, new_column_names[1:length(larger_values)]] <- larger_values
  } else {
    # If GaussianDOYLate is not a finite number, fill the corresponding columns with NA
    SpeciesSlope[i, new_column_names] <- NA
  }
}

# Check the updated SpeciesSlope data frame
head(SpeciesSlope)


#Pivot the table (new values in single column)

# Pivot the data frame
SpeciesSlope_pivoted <- pivot_longer(
  data = SpeciesSlope,
  cols = 19:38,
  names_to = "Position",
  values_to = "Difference"
)

# Merge SpeciesSlope_pivoted with GaussianCurve based on matching values in Difference and x
merged_data <- merge(SpeciesSlope_pivoted, GaussianCurve, by.x = "Difference", by.y = "x", all.x = TRUE)

# Fill in RelativeAbundance with density values where available
merged_data$RelativeAbundance[!is.na(merged_data$density)] <- merged_data$density[!is.na(merged_data$density)]

# Drop the density column as it's no longer needed
merged_data <- merged_data[, !(names(merged_data) %in% "density")]


# Add the DOY values for the data
# Initialize DOY column
merged_data$DOY <- NA

# Fill in DOY values based on Position
merged_data$DOY <- ifelse(
  substr(merged_data$Position, 1, 4) == "DOYE",
  merged_data$DOY1 - as.numeric(substr(merged_data$Position, 5, 6)),
  merged_data$DOY4 + as.numeric(substr(merged_data$Position, 5, 6))
)

# Add values on the difference between ranges of photos and flower observations
FloweringEndPivot
FloweringStartPivot


# Join to merge FloweringPhotoStart into merged_data
merged_data <- merged_data %>%
  left_join(FloweringStartPivot %>% select(camera, Species, FloweringPhotoStart),
            by = c("camera", "Species"))

# Join to merge FloweringPhotoEnd into merged_data
merged_data <- merged_data %>%
  left_join(FloweringEndPivot %>% select(camera, Species, FloweringPhotoEnd),
            by = c("camera", "Species"))


# Calculate the flower abundances
# Conditionally standardised to be 1 on the DOY that extrapolation starts
merged_data$RelativeAbundanceX <- with(merged_data,
                                       ifelse(
                                         substr(Position, 1, 4) %in% c("DOYE", "DOYL"),
                                         RelativeAbundance /
                                           ave(RelativeAbundance,
                                               camera,
                                               Species,
                                               substr(Position, 1, 4),
                                               FUN = function(x) max(x, na.rm = TRUE)
                                           ),
                                         NA_real_
                                       )
)

# Calculate the predicted abundances
merged_data$PredictedAbundance <- with(merged_data, ifelse(
  grepl("^DOYE", Position),
  RelativeAbundanceX * AbundanceEarly,
  ifelse(grepl("^DOYL", Position),
         RelativeAbundanceX * AbundanceLate,
         NA  # You can set a default value if needed
  )
))


# Save the file
write.csv(merged_data,'./data/GaussianPredictions.csv')

# Merge the Gaussian extrapolations to CrosstabMax
# Load and edit GaussianPredictions for join
GaussianPredictions <- read.csv('./data/GaussianPredictions.csv')
#Remove "X" column
GaussianPredictions <- GaussianPredictions[, !(names(GaussianPredictions) %in% c('X'))]

# Replacing values with NA in PredictedAbundance based on conditions
GaussianPredictions$PredictedAbundance <- ifelse(
  GaussianPredictions$FloweringPhotoStart > 3 & grepl("^DOYE", GaussianPredictions$Position),
  NA,
  GaussianPredictions$PredictedAbundance
)

# Replacing values with NA in PredictedAbundance based on conditions
GaussianPredictions$PredictedAbundance <- ifelse(
  GaussianPredictions$FloweringPhotoEnd > 3 & grepl("^DOYL", GaussianPredictions$Position),
  NA,
  GaussianPredictions$PredictedAbundance
)

# Removing rows where PredictedAbundance is NA
GaussianPredictions <- GaussianPredictions[!is.na(GaussianPredictions$PredictedAbundance), ]


# Load and edit CrosstabMax for join
CrosstabMax <- read.csv('./data/CrosstabMax.csv')
#Remove "X" column
CrosstabMax <- CrosstabMax[, !(names(CrosstabMax) %in% c('X'))]

# Pivot the table for join
CrosstabMaxPivoted <- CrosstabMax %>%
  pivot_longer(cols = 4:36, names_to = "Species", values_to = "Abundance")

# Prepare GaussianPredictions for merging by selecting the relevant columns and renaming them
GaussianPredictionsPrepared <- GaussianPredictions %>%
  select(camera, DOY, Species, Abundance = PredictedAbundance)

# Combine the two data frames
CrosstabMaxCombined <- bind_rows(CrosstabMaxPivoted, GaussianPredictionsPrepared)

# Save the file
write.csv(CrosstabMaxCombined,'./data/FloweringANDExtrapolatedData.csv')

# Then the interpolation of the values based on the observed values and predicted values
# Load and edit
FloweringANDExtrapolatedData <- read.csv('./data/FloweringANDExtrapolatedData.csv')
#Remove "X" column
FloweringANDExtrapolatedData <- FloweringANDExtrapolatedData[, !(names(FloweringANDExtrapolatedData) %in% c('X'))]

# Extend the dataframe so that it has all the DOY values for all the cameras
# Create a data frame with all combinations of camera, Species, and DOY (+- of the DOY range of FloweringANDExtrapolatedData)
all_combinations <- expand.grid(
  camera = unique(FloweringANDExtrapolatedData$camera),
  Species = unique(FloweringANDExtrapolatedData$Species),
  DOY = 87:289
)

# Merge the original data frame with the expanded data frame
FloweringANDExtrapolatedData_extended <- merge(all_combinations, FloweringANDExtrapolatedData,
                                               by = c("camera", "Species", "DOY"), all.x = TRUE)

# Replace all zeros with NA in Abundance column
FloweringANDExtrapolatedData_extended$Abundance[FloweringANDExtrapolatedData_extended$Abundance == 0] <- NA

# Add a column of 0 for calculating moving averages
# Group by camera and Species, calculate the minimum non-zero DOY and maximum non-zero DOY, and create the Calculation0 column

FloweringANDExtrapolatedData_extended <- FloweringANDExtrapolatedData_extended %>%
  group_by(camera, Species) %>%
  mutate(
    min_nonzero_DOY = ifelse(any(Abundance != 0, na.rm = TRUE), min(DOY[Abundance != 0], na.rm = TRUE), NA),
    max_nonzero_DOY = ifelse(any(Abundance != 0, na.rm = TRUE), max(DOY[Abundance != 0], na.rm = TRUE), NA),
    StartEnd0 = ifelse((is.na(min_nonzero_DOY) | DOY < min_nonzero_DOY) | (is.na(max_nonzero_DOY) | DOY > max_nonzero_DOY), 0, NA)
  ) %>%
  select(-min_nonzero_DOY, -max_nonzero_DOY)  # Remove the intermediate columns


# Add 0 for the gaps
# Function to calculate GapSize based on identified gaps
calculate_gap_size <- function(abundance, doy) {
  gap_sizes <- rep(0, length(abundance))
  start_index <- NA
  for (i in seq_along(abundance)) {
    if (!is.na(abundance[i]) && abundance[i] != 0) {
      if (!is.na(start_index)) {
        # Skip writing down the first value of each gap
        gap_sizes[(start_index+1):(i-1)] <- i - start_index - 1
      }
      start_index <- i
    }
  }
  return(gap_sizes)
}

# Group by camera and Species, calculate GapSize for each group
FloweringANDExtrapolatedData_extended <- FloweringANDExtrapolatedData_extended %>%
  group_by(camera, Species) %>%
  mutate(
    GapSize = calculate_gap_size(Abundance, DOY)
  )

# Write 0 to Gap0 if GapSize is 5 or larger
FloweringANDExtrapolatedData_extended <- FloweringANDExtrapolatedData_extended %>%
  mutate(
    Gap0 = ifelse(GapSize >= 5, 0, NA)
  )

# Combine values from Gap0 and StartEnd0 into Calculation0 by taking the maximum value
FloweringANDExtrapolatedData_extended <- FloweringANDExtrapolatedData_extended %>%
  mutate(
    Calculation0 = pmax(Gap0, StartEnd0, na.rm = TRUE)
  )


# Function to perform linear interpolation
linear_interpolate <- function(abundance, doy) {
  na_values <- is.na(abundance)
  if (any(!na_values)) {
    abundance <- na.approx(abundance, x = doy, na.rm = FALSE)
  }
  return(abundance)
}

# Calculate AbundanceInterpolation based on the given conditions
FloweringANDExtrapolatedData_extended <- FloweringANDExtrapolatedData_extended %>%
  mutate(
    AbundanceInterpolation = ifelse(Abundance != 0, Abundance, NA),  # Step 1
    AbundanceInterpolation = ifelse(GapSize < 5, linear_interpolate(AbundanceInterpolation, DOY), AbundanceInterpolation)  # Step 2
  )

# Combine values from Calculation0 and AbundanceInterpolation into AbundanceInterpolationCalculation0 by taking the maximum value
FloweringANDExtrapolatedData_extended <- FloweringANDExtrapolatedData_extended %>%
  mutate(
    AbundanceInterpolationCalculation0 = pmax(Calculation0, AbundanceInterpolation, na.rm = TRUE)
  )

# Interpolate non-zero values to 0 in larger gaps in data
# Function to identify adjacent 0 and non-zero values and perform interpolation
interpolate_gaps <- function(doy, abundance) {
  gap_interpolate <- rep(NA, length(abundance))
  n <- length(abundance)

  for (i in 1:(n - 1)) {
    if (abundance[i] == 0 && abundance[i + 1] != 0) {
      gap_interpolate[i] <- abundance[i + 1] * 2 / 3
      if (i > 1 && abundance[i - 1] == 0) {
        gap_interpolate[i - 1] <- abundance[i + 1] * 1 / 3
      }
    } else if (abundance[i] != 0 && abundance[i + 1] == 0) {
      gap_interpolate[i + 1] <- abundance[i] * 2 / 3
      if (i < n - 1 && abundance[i + 2] == 0) {
        gap_interpolate[i + 2] <- abundance[i] * 1 / 3
      }
    }
  }

  return(gap_interpolate)
}

# Apply the function to the data frame for each combination of camera and Species
FloweringANDExtrapolatedData_extended <- FloweringANDExtrapolatedData_extended %>%
  group_by(camera, Species) %>%
  mutate(
    GapInterpolate = interpolate_gaps(DOY, AbundanceInterpolationCalculation0)
  ) %>%
  ungroup()

# Calculate sum of values for each row
FloweringANDExtrapolatedData_extended <- FloweringANDExtrapolatedData_extended %>%
  mutate(
    PredictedAbundance = rowSums(select(., AbundanceInterpolationCalculation0, GapInterpolate), na.rm = TRUE)
  )

# Save the file
write.csv(FloweringANDExtrapolatedData_extended,'./data/PredictedFlowerAbundance.csv')

#Create a Crosstab out of the predicted flower abundances
PredictedFlowerAbundance <- read.csv('./data/PredictedFlowerAbundance.csv')
#Remove unnecessary columns
PredictedFlowerAbundance <- PredictedFlowerAbundance[, !(names(PredictedFlowerAbundance) %in% c('X'))]
PredictedFlowerAbundance <- PredictedFlowerAbundance[, !(names(PredictedFlowerAbundance) %in% c('filedate'))]
PredictedFlowerAbundance <- PredictedFlowerAbundance[, !(names(PredictedFlowerAbundance) %in% c('StartEnd0'))]
PredictedFlowerAbundance <- PredictedFlowerAbundance[, !(names(PredictedFlowerAbundance) %in% c('GapSize'))]
PredictedFlowerAbundance <- PredictedFlowerAbundance[, !(names(PredictedFlowerAbundance) %in% c('Gap0'))]
PredictedFlowerAbundance <- PredictedFlowerAbundance[, !(names(PredictedFlowerAbundance) %in% c('ele_m'))]


# Create the crosstab
# Summarise duplicates of PredictedAbundance
PredictedFlowerAbundance <- PredictedFlowerAbundance %>%
  group_by(camera, DOY, Species) %>%
  summarise(PredictedAbundance = sum(PredictedAbundance, na.rm = TRUE), .groups = 'drop')

# Create the crosstab
PredictedFlowerAbundanceCrosstab <- PredictedFlowerAbundance %>%
  pivot_wider(
    names_from = Species,
    values_from = PredictedAbundance,
    values_fill = list(PredictedAbundance = 0) # Fill missing values with 0
  )

# Add the meta data
# Camera metadata
reconmeta <- read_csv('./data/reconycx_trailcam_metadata.csv')
reconmetadata  <- reconmeta %>% mutate(camera= paste0("cam",camera) ) %>% as.data.frame()

# Merge the reconmetadata columns to PredictedFlowerAbundanceCrosstab
PredictedFlowerAbundanceCrosstab <- merge(PredictedFlowerAbundanceCrosstab, reconmetadata[, c("camera", "habitat", "elev_m", "canopy_cov")], by = "camera")

# Save the file
write.csv(PredictedFlowerAbundanceCrosstab,'./data/PredictedFlowerAbundanceCrosstab.csv')



# Calculating the earliest, latest, and flowering season length for all the species at all the cameras####

#Formatting the prediction data frame (this is created in the last part of the script)
CrosstabMax <- read.csv('./data/PredictedFlowerAbundanceCrosstab.csv')
names(CrosstabMax)[names(CrosstabMax) == "elev_m.x"] <- "elev_m"
CrosstabMax <- CrosstabMax[, !(names(CrosstabMax) %in% c('habitat.x'))]
CrosstabMax <- CrosstabMax[, !(names(CrosstabMax) %in% c('habitat.y'))]
CrosstabMax <- CrosstabMax[, !(names(CrosstabMax) %in% c('elev_m.y'))]
CrosstabMax <- CrosstabMax[, !(names(CrosstabMax) %in% c('canopy_cov'))]
CrosstabMax <- CrosstabMax[, c(names(CrosstabMax)[1:3], "elev_m", names(CrosstabMax)[4:(ncol(CrosstabMax)-1)])]

head(CrosstabMax)

#Remove "X" column
CrosstabMax <- CrosstabMax[, !(names(CrosstabMax) %in% c('X'))]

# Only keeping the earliest DOY for all the species
species_columns <- colnames(CrosstabMax[4:36])

# Create an empty data frame to store the results
CrosstabEarly <- unique(CrosstabMax[, c("camera", "elev_m")])

# Iterate over each species
for (species in species_columns) {
  # Find the smallest DOY for the current species for each camera and elev_m combination
  min_doy_species <- aggregate(DOY ~ camera + elev_m, data = subset(CrosstabMax, CrosstabMax[[species]] > 0), FUN = min)

  # Merge the result with the original data frame using a left join
  CrosstabEarly <- merge(CrosstabEarly, min_doy_species, by.x = c("camera", "elev_m"), by.y = c("camera", "elev_m"), all.x = TRUE)

  # Rename the new column with the original species name
  colnames(CrosstabEarly)[ncol(CrosstabEarly)] <- species
}


# Only keeping the latest DOY for all the species

# Create an empty data frame to store the results
CrosstabLate <- unique(CrosstabMax[, c("camera", "elev_m")])

# Iterate over each species
for (species in species_columns) {
  # Find the smallest DOY for the current species for each camera and elev_m combination
  max_doy_species <- aggregate(DOY ~ camera + elev_m, data = subset(CrosstabMax, CrosstabMax[[species]] > 0), FUN = max)

  # Merge the result with the original data frame using a left join
  CrosstabLate <- merge(CrosstabLate, max_doy_species, by.x = c("camera", "elev_m"), by.y = c("camera", "elev_m"), all.x = TRUE)

  # Rename the new column with the original species name
  colnames(CrosstabLate)[ncol(CrosstabLate)] <- species
}


#Calculating the flowering season lenght for all the species
subtracted_values <- CrosstabLate[, 3:35] - CrosstabEarly[, 3:35]

# Create a new data frame CrosstabSeasonLength
CrosstabSeasonLength <- data.frame(
  CrosstabLate$camera,
  CrosstabLate$elev_m,
  subtracted_values
)

# Calculating and plotting the average flowering phenology for the species -------------

#Calculating the total number of flowers
CrosstabMax <- read.csv('./data/PredictedFlowerAbundanceCrosstab.csv')
head(CrosstabMax)
#Remove "X" column
CrosstabMax <- CrosstabMax[, !(names(CrosstabMax) %in% c('X'))]

# Exclude unnecessary columns before grouping and summarizing
columns_to_remove <- c('DOY', 'habitat', 'canopy_cov')
CrosstabMax_filtered <- CrosstabMax[, !(names(CrosstabMax) %in% columns_to_remove)]

# Summing things
NumberOfFlowers <- CrosstabMax_filtered %>%
  group_by(camera, elev_m) %>%
  summarise_all(sum)


# Calculating the sums of flower number * DOY
# Exclude unnecessary columns before grouping and summarizing
columns_to_remove <- c('habitat', 'canopy_cov')
CrosstabMax_filtered <- CrosstabMax[, !(names(CrosstabMax) %in% columns_to_remove)]
CrosstabMax_filtered <- CrosstabMax_filtered[, c(setdiff(names(CrosstabMax_filtered), "DOY"), "DOY")]

# Multiply values of selected columns by the 'DOY' column
for (i in 2:34) {
  CrosstabMax_filtered[, i] <- CrosstabMax_filtered[, i] * CrosstabMax_filtered[, 36]
}

#Summing things
SumFlowerDOY <- CrosstabMax_filtered %>%
  group_by(camera, elev_m) %>%
  summarise_all(sum)

#Remove DOY column
SumFlowerDOY <- SumFlowerDOY[1:35]

#Calculating the average flowering dates (sum(flower number * DOY)/total number of flowers)
average_flowering_values <- SumFlowerDOY[, 3:35] / NumberOfFlowers[, 3:35]

# Create a new data frame AverageFlowering
AverageFlowering <- data.frame(
  CrosstabLate$camera,
  CrosstabLate$elev_m,
  average_flowering_values
)


# Counting local species number and flower count --------------------------
#Counting species number
NumberOfFlowers
columns_to_remove <- c('elev_m', 'X')
NumberOfFlowers <- NumberOfFlowers[, !(names(NumberOfFlowers) %in% columns_to_remove)]

# Exclude unnecessary columns
selected_columns <- NumberOfFlowers[,2:34]

# Calculate the species number
SpeciesNumber <- data.frame(
  camera = NumberOfFlowers$camera,
  SpeciesNumber = rowSums(selected_columns != 0))

# Calculate the flower number
FlowerNumber <- data.frame(
  camera = NumberOfFlowers$camera,
  FlowerNumber = rowSums(selected_columns))


# Collect interesting numbers from the previous parts -------

#Local averages of earliest and latest flowering dates and flowering season lengths of all species

# THIS FOLLOWING IS ONLY TO CHECK THE SENSITIVITY OF OBSERVATION PERIOD ON RESULTS, FOR MAIN RESULTS SKIP
# THIS SECTION FILTERS OUT THE CASES WHEN FLOWERING MEETS THE EDGE OF THE OBSERVATION PERIOd
# Identifying species hitting the observation period borders####
FloweringStart
TooEarly <- FloweringStart[ , -2]
TooEarly[ , -1] <- ifelse(TooEarly[ , -1] < 0, 1, 0)

write.csv(TooEarly,'./data/TooEarlyFilter.csv')

FloweringEnd
TooLate <- FloweringEnd[ , -2]
TooLate[ , -1] <- ifelse(TooLate[ , -1] < 0, 1, 0)

write.csv(TooLate,'./data/TooLateFilter.csv')

TooEarly <- read.csv('./data/TooEarlyFilter.csv')
TooLate <- read.csv('./data/TooLateFilter.csv')

# Remove earliest dates
# 1. Reorder TooEarly to match the row order of CrosstabEarly by 'camera'
TooEarly_ord <- TooEarly[ match(CrosstabEarly$camera, TooEarly$camera), ]

# 2. Identify which columns to process (all but 'camera')
common_cols <- setdiff(intersect(names(CrosstabEarly), names(TooEarly)), "camera")

# 3. Loop over each common column and set CrosstabEarly to NA where TooEarly == 1
for (col in common_cols) {
  # logical index of TooEarly == 1
  to_na <- TooEarly_ord[[col]] == 1

  # assign NA in CrosstabEarly
  CrosstabEarly[[col]][ to_na ] <- NA
}

# Remove the latest dates
TooLate_ord <- TooLate[ match(CrosstabLate$camera, TooLate$camera), ]

# 2. Identify which columns to process (all but 'camera')
common_cols <- setdiff(intersect(names(CrosstabLate), names(TooLate)), "camera")

# 3. Loop over each common column and set CrosstabLate to NA where TooLate == 1
for (col in common_cols) {
  # logical index of TooEarly == 1
  to_na <- TooLate_ord[[col]] == 1

  # assign NA in CrosstabEarly
  CrosstabLate[[col]][ to_na ] <- NA
}

# Remove season length
# Because of earliest date
# 1. Reorder TooEarly to match the row order of CrosstabEarly by 'camera'
TooEarlySeasonLength_ord <- TooEarly[ match(CrosstabSeasonLength$camera, TooEarly$camera), ]

# 2. Identify which columns to process (all but 'camera')
common_cols <- setdiff(intersect(names(CrosstabSeasonLength), names(TooEarly)), "camera")

# 3. Loop over each common column and set CrosstabEarly to NA where TooEarly == 1
for (col in common_cols) {
  # logical index of TooEarly == 1
  to_na <- TooEarlySeasonLength_ord[[col]] == 1

  # assign NA in CrosstabEarly
  CrosstabSeasonLength[[col]][ to_na ] <- NA
}

# Because of latest date
# 1. Reorder TooEarly to match the row order of CrosstabEarly by 'camera'
TooLateSeasonLength_ord <- TooLate[ match(CrosstabSeasonLength$camera, TooLate$camera), ]

# 2. Identify which columns to process (all but 'camera')
common_cols <- setdiff(intersect(names(CrosstabSeasonLength), names(TooLate)), "camera")

# 3. Loop over each common column and set CrosstabEarly to NA where TooEarly == 1
for (col in common_cols) {
  # logical index of TooEarly == 1
  to_na <- TooLateSeasonLength_ord[[col]] == 1

  # assign NA in CrosstabEarly
  CrosstabSeasonLength[[col]][ to_na ] <- NA
}

# THE REMOVAL PART ENDS HERE

AverageEarly <- data.frame(camera = CrosstabEarly[, 1],
                           Early = rowMeans(CrosstabEarly[, 3:35], na.rm = TRUE))
AverageLate <- data.frame(camera = CrosstabLate[, 1],
                          Late = rowMeans(CrosstabLate[, 3:35], na.rm = TRUE))
AverageAverageFlowering <- data.frame(camera = AverageFlowering[, 1],
                          AverageFlowering = rowMeans(AverageFlowering[, 3:35], na.rm = TRUE))
AverageSeasonLength <- data.frame(camera = CrosstabSeasonLength[, 1],
                                  SeasonLength = rowMeans(CrosstabSeasonLength[, 3:35], na.rm = TRUE))

EarlyLateSeason <- reconmetadata %>%
  left_join(SpeciesNumber,
            by = "camera") %>%
  left_join(FlowerNumber,
            by = "camera") %>%
  left_join(AverageEarly,
            by = "camera") %>%
  left_join(AverageLate,
            by = "camera") %>%
  left_join(AverageAverageFlowering,
            by = "camera") %>%
  left_join(AverageSeasonLength,
            by = "camera") %>%
    left_join(NumberOfFlowers,
            by = "camera")

write.csv(EarlyLateSeason,'./data/EarlyLateSeason.csv')

# Species level patterns --------------------------------------------------

#Best species (abundant, in many cameras)
#1. Ranunculus_montanus (all sites, 9 sites with less than 10)
#2. Crocus_albiflorus (all sites, 9 sites with less than 10)
#3. Alchemilla_conjuncta (all sites, 12 sites with less than 10)
#4. Campanula_scheuchzeri (34/35 sites, 15 sites with less than 10)
#5. Anthyllis_vulneraria (33/35 sites, 18 sites with less than 10)
#6. Alchemilla_xanthochlora (30/35 sites, 22 sites with less than 10)

# Checking all the pieces
AverageFlowering
colnames(AverageFlowering)[colnames(AverageFlowering) == "CrosstabLate.camera"] <- "camera"
colnames(AverageFlowering)[colnames(AverageFlowering) == "CrosstabLate.elev_m"] <- "elev_m"

CrosstabEarly
CrosstabLate
CrosstabSeasonLength
colnames(CrosstabSeasonLength)[colnames(CrosstabSeasonLength) == "CrosstabLate.camera"] <- "camera"
colnames(CrosstabSeasonLength)[colnames(CrosstabSeasonLength) == "CrosstabLate.elev_m"] <- "elev_m"

# Combining and saving the data
SpeciesTraits <- reconmetadata %>%
  left_join(AverageFlowering,
            by = "camera") %>%
  left_join(CrosstabEarly,
            by = "camera") %>%
  left_join(CrosstabLate,
            by = "camera") %>%
  left_join(CrosstabSeasonLength,
            by = "camera")


# Change "alpine" habitat to "open"
SpeciesTraits$habitat[SpeciesTraits$habitat == "alpine"] <- "open"

SpeciesTraits$site = as.factor(SpeciesTraits$site)
levels(SpeciesTraits$site)[levels(SpeciesTraits$site) == "nbo"] <- "below treeline"
levels(SpeciesTraits$site)[levels(SpeciesTraits$site) == "subcal"] <- "at treeline"
levels(SpeciesTraits$site)[levels(SpeciesTraits$site) == "cal"] <- "above treeline"
levels(SpeciesTraits$site)

write.csv(SpeciesTraits,'./data/SpeciesTraits.csv')
#Note: mean x, first y, last x.x, season length y.y


#Making the GLMs for phenology mean, first, last, season length
SpeciesTraits <- read.csv('./data/SpeciesTraits.csv')

#GLMs for mean flowering
columns_to_analyze <- c(
  "Ranunculus_montanus.x",
  "Crocus_albiflorus.x",
  "Alchemilla_conjuncta.x",
  "Campanula_scheuchzeri.x",
  "Anthyllis_vulneraria.x",
  "Alchemilla_xanthochlora.x"
)

for (i in 1:length(columns_to_analyze)) {
  formula_str <- paste(columns_to_analyze[i], "~ site * habitat")
  glm_name <- paste("M3.", i, sep = "")

  assign(glm_name, glm(formula_str, data = SpeciesTraits))

  cat("GLM", glm_name, "created and run:\n")
  print(get(glm_name))
  cat("\n")
}

#GLMs for earliest flowering
for (i in 1:length(columns_to_analyze)) {
  formula_str <- paste(columns_to_analyze[i], "~ site * habitat")
  glm_name <- paste("M4.", i, sep = "")

  assign(glm_name, glm(formula_str, data = SpeciesTraits))

  cat("GLM", glm_name, "created and run:\n")
  print(get(glm_name))
  cat("\n")
}


#GLMs for latest flowering
for (i in 1:length(columns_to_analyze)) {
  formula_str <- paste(columns_to_analyze[i], "~ site * habitat")
  glm_name <- paste("M5.", i, sep = "")

  assign(glm_name, glm(formula_str, data = SpeciesTraits))

  cat("GLM", glm_name, "created and run:\n")
  print(get(glm_name))
  cat("\n")
}


#GLMs for season length
for (i in 1:length(columns_to_analyze)) {
  formula_str <- paste(columns_to_analyze[i], "~ site * habitat")
  glm_name <- paste("M6.", i, sep = "")

  assign(glm_name, glm(formula_str, data = SpeciesTraits))

  cat("GLM", glm_name, "created and run:\n")
  print(get(glm_name))
  cat("\n")
}


# Making the phenology plots
# Define the desired order of levels
desired_order <- c("below treeline.open", "below treeline.canopy", "at treeline.open", "at treeline.canopy","above treeline.open")

# Map the levels to custom labels
custom_labels <- c("below treeline, open", "below treeline, canopy", "at treeline, open", "at treeline, canopy","above treeline")

# Create a new factor with the reordered levels
SpeciesTraits$interaction <- factor(interaction(SpeciesTraits$site, SpeciesTraits$habitat))

# Average flowering phenology
# Define the species columns to include
species_columns <- c("Ranunculus_montanus.x", "Crocus_albiflorus.x", "Alchemilla_conjuncta.x",
                     "Campanula_scheuchzeri.x", "Anthyllis_vulneraria.x", "Alchemilla_xanthochlora.x")

# Reshape the data for faceting
melted_data <- melt(SpeciesTraits, id.vars = c("interaction"), measure.vars = species_columns)

# If 'interaction' is not already a factor, or if you want to set a specific order:
desired_order <- c("below treeline.open", "below treeline.canopy", "at treeline.open", "at treeline.canopy","above treeline.open")  # Replace with your desired order
melted_data$interaction <- factor(melted_data$interaction, levels = desired_order)

# Define custom labels for the x-axis
custom_labels <- c("below treeline, open", "below treeline, canopy", "at treeline, open", "at treeline, canopy","above treeline")  # Replace with your custom labels

# Define custom legend title and labels
legend_title <- "Plant Species"
legend_title_flip <- "Elevation & Habitat"
custom_legend_labels <- c("Ranunculus montanus", "Crocus albiflorus", "Alchemilla conjuncta", "Campanula scheuchzeri", "Anthyllis vulneraria", "Alchemilla xanthochlora")

# Create the plot
####Supplement FigDd####
SFigDd <- ggplot(melted_data, aes(x = interaction, y = value, fill = variable)) +
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.7) +
  #geom_jitter(position = position_dodge(width = 0.8), alpha = 0.6, size = 1) +
  labs(title = "Peak Flowering",
       x = "Elevation",
       y = "Day of Year",,
       fill = legend_title) +
  scale_x_discrete(labels = custom_labels) +
  scale_fill_manual(values = c("yellow2", "thistle3", "yellowgreen", "slateblue3", "goldenrod1", "springgreen"),
                    labels = custom_legend_labels) + # Set custom labels
  coord_flip() +  # Flip x and y axes
  theme_minimal(base_size = 24)

SFigDd

# Earliest flowering phenology
# Define the species columns to include
species_columns <- c("Ranunculus_montanus.y", "Crocus_albiflorus.y", "Alchemilla_conjuncta.y",
                     "Campanula_scheuchzeri.y", "Anthyllis_vulneraria.y", "Alchemilla_xanthochlora.y")

# Reshape the data for faceting
melted_data <- melt(SpeciesTraits, id.vars = c("interaction"), measure.vars = species_columns)

####Supplement FigDa####
# Plot the original values of specified species using a boxplot with data points and facet by interaction
SFigDa <- ggplot(melted_data, aes(x = interaction, y = value, fill = variable)) +
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.7) +
  #geom_jitter(position = position_dodge(width = 0.8), alpha = 0.6, size = 1) +
  labs(title = "First Flowering",
       x = "Elevation",
       y = "Day of Year",,
       fill = legend_title) +
  scale_x_discrete(labels = custom_labels) +
  scale_fill_manual(values = c("yellow2", "thistle3", "yellowgreen", "slateblue3", "goldenrod1", "springgreen"),
                    labels = custom_legend_labels) + # Set custom labels
  coord_flip() +  # Flip x and y axes
  theme_minimal(base_size = 24)
SFigDa

# Latest flowering phenology
# Define the species columns to include
species_columns <- c("Ranunculus_montanus.x.x", "Crocus_albiflorus.x.x", "Alchemilla_conjuncta.x.x",
                     "Campanula_scheuchzeri.x.x", "Anthyllis_vulneraria.x.x", "Alchemilla_xanthochlora.x.x")

# Reshape the data for faceting
melted_data <- melt(SpeciesTraits, id.vars = c("interaction"), measure.vars = species_columns)

####Supplement FigDb####
# Plot the original values of specified species using a boxplot with data points and facet by interaction
SFigDb <- ggplot(melted_data, aes(x = interaction, y = value, fill = variable)) +
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.7) +
  #geom_jitter(position = position_dodge(width = 0.8), alpha = 0.6, size = 1) +
  labs(title = "Last Flowering",
       x = "Elevation",
       y = "Day of Year",,
       fill = legend_title) +
  scale_x_discrete(labels = custom_labels) +
  scale_fill_manual(values = c("yellow2", "thistle3", "yellowgreen", "slateblue3", "goldenrod1", "springgreen"),
                    labels = custom_legend_labels) + # Set custom labels
  coord_flip() +  # Flip x and y axes
  theme_minimal(base_size = 24)
SFigDb

# Season length
# Define the species columns to include
species_columns <- c("Ranunculus_montanus.y.y", "Crocus_albiflorus.y.y", "Alchemilla_conjuncta.y.y",
                     "Campanula_scheuchzeri.y.y", "Anthyllis_vulneraria.y.y", "Alchemilla_xanthochlora.y.y")

# Reshape the data for faceting
melted_data <- melt(SpeciesTraits, id.vars = c("interaction"), measure.vars = species_columns)

####Supplement FigDc####
# Plot the original values of specified species using a boxplot with data points and facet by interaction
SFigDc <-ggplot(melted_data, aes(x = interaction, y = value, fill = variable)) +
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.7) +
  #geom_jitter(position = position_dodge(width = 0.8), alpha = 0.6, size = 1) +
  labs(title = "Season Length",
       x = "Elevation",
       y = "Flowering Duration (Days)",,
       fill = legend_title) +
  scale_x_discrete(labels = custom_labels) +
  scale_fill_manual(values = c("yellow2", "thistle3", "yellowgreen", "slateblue3", "goldenrod1", "springgreen"),
                    labels = custom_legend_labels) + # Set custom labels
  coord_flip() +  # Flip x and y axes
  theme_minimal(base_size = 24)
SFigDc

# Combine the plots into a single layout
combined_plot_FigSD <- (SFigDa | SFigDb) / (SFigDc | SFigDd) + labs(title = "") +plot_layout(guides = "collect") & plot_annotation(tag_levels = "a")  &
  theme(legend.position = "bottom")

# Display the combined plot
print(combined_plot_FigSD)

ggsave("./figs/SupplementFigureD.png", width = 40, height = 40, units = "cm")
