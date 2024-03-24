rm(list=ls())

library(readxl)
library(dplyr)
# read AADR dataset
AADR <- read_excel("/Users/tianmi/Desktop/Bioinf/BINP29/Population/AADR/AADR_Annotation.xlsx")
# View(AADR)

# Read a list of IDs from Y haplogroups from a text file
ID_list <- read.table("/Users/tianmi/Desktop/Bioinf/BINP29/Population/AncientYDNA/ID_once.txt", header = FALSE)

# Read the Y_ISOGG file, which contains information about the relationship of Y haplogroups
Y_ISOGG <- read.delim("/Users/tianmi/Desktop/Bioinf/BINP29/Population/AncientYDNA/chrY_hGrpTree_isogg2016.txt", header = FALSE)

# Extract relevant columns of Y haplogroup data from the AADR dataset
Y_data <- AADR[,c(2, 3, 9, 15, 27)]

# Check and revise column names for clarity
names(Y_data)
names(Y_data) <- c("Genetic_ID", "Master_ID", "Full_Date", "Political_Entity", "Y_haplo")
nrow(Y_data) # 16466
Y_data$Political_Entity <- gsub("Gernamy", "Germany", Y_data$Political_Entity)

Y_data <- Y_data %>%
  mutate(
    Y_haplo = gsub("\\(.*?\\)", "", Y_haplo),
    Y_haplo = gsub("[-_~*;x?! ].*", "", Y_haplo),
    Y_haplo = gsub("^\\s+|\\s+$", "", Y_haplo)
  )

# final cleaned data frame 
cleaned_Y <- Y_data %>%
  filter(grepl("^[A-T]", Y_haplo)) %>%
  select(-Genetic_ID, -Master_ID) %>%
  select("Political_Entity", "Y_haplo", "Full_Date")
  

### Calculate the frequency of each Y haplogroup within each country
country_freq <- Y_data %>% 
  group_by(Political_Entity, Y_haplo) %>%
  summarize(country_count = n(), .groups = 'drop') %>% 
  group_by(Political_Entity) %>%
  # Calculate the frequency of haplogroups within each country
  mutate(country_frequency = country_count / sum(country_count)) %>%
  ungroup()

# Delete and correct misformatted haplogroup ID 
country_freq <- country_freq %>%
  filter(grepl("^[A-T]", Y_haplo)) 

### Calculate the frequency of each Y haplogroup in the world
global_freq <- Y_data %>%
  group_by(Y_haplo) %>%
  summarise(global_count = n(), .groups = 'drop') %>%
  mutate(global_frequency = global_count / sum(global_count))

# Delete and correct misformatted haplogroup ID 
global_freq <- global_freq %>%
  filter(grepl("^[A-T]", Y_haplo)) 

### Categorize countries
# Categorize data based on countries 
# Calculate how many different haplogroups there are in each country and how many haplogroups there are in total
country_hap <- country_freq %>% group_by(Political_Entity) %>% 
  summarize(type_haps = n_distinct(Y_haplo),
            total_num_haps = sum(country_count))

options(scipen = 999)
Population_vector <- c(45810000, 2800000, 8950000, 10140000, 408000,
                       11600000, 12000000, 2600000, 27200000, 2200000,
                       3703, 1400000000, 4000000, 11260000, 10510000,
                       96000000, 5800000, 1300000, 60000, 65000000,
                       83000000, 60000, 400000, 10000000, 400000,
                       88000000, 5100000, 100000, 9500000, 60000000,
                       130000000, 11200000, 19000000, 54000000, 6700000,
                       2000000, 20000000, 3500000, 30050000, 18000000,
                       5500000, 34000000, 38000000, 10500000, 3400000,
                       20000000, 150000000, 60000000, 50000000, 200000,
                       11000000, 8800000, 22000000, 24000000, 65000000,
                       85000000, 6500000, 350000000, 45000000, 68000000)
country_hap$Population <- Population_vector

# Merge country info
merged_data <- merge(cleaned_Y, country_freq, by=c("Political_Entity", "Y_haplo"), all = TRUE) 
merged_data <- merge(merged_data, country_hap, by="Political_Entity", all = TRUE) %>%
  rename(country_population = Population)
merged_data <- merged_data %>% 
  # Estimate the number of individuals for each haplogroup by multiplying the frequency by population, then rounding to the nearest whole number
  mutate(country_individuals = round(country_frequency * country_population, digits = 0))

# Merge global info
merged_data <- merge(merged_data, global_freq, by="Y_haplo", all = TRUE) 
merged_data$global_population <- 7888000000
merged_data <- merged_data %>% 
  mutate(global_individuals = round(global_frequency * global_population, digits = 0))

# Create a file with unique haplogroup IDs
unique_haplogroups <- unique(cleaned_Y$Y_haplo)
unique_haplogroups_df <- data.frame(Unique_Haplogroups = unique_haplogroups)
write.csv(unique_haplogroups_df, "unique_haplogroups.csv", row.names = FALSE)

### Haplogroup Path
# Create a df to store relationships between haplogroups (father-son relationships)
relationships <- data.frame(
  father = Y_ISOGG$V2, # Assuming column V2 contains father haplogroups
  son = Y_ISOGG$V1 # Assuming column V1 contains son haplogroups
)

### Find haplogroup path
find_path <- function(current, relationships) {
  # # Initialize the path with the current haplogroup
  path <- current
  # Find the father of the current haplogroup where current haplogroup is the son and ensure the father is not a placeholder (e.g., "#")
  father <- relationships$father[relationships$son == current & relationships$father !="#"]
  # # If a father exists, recursively call find_path to construct the path from the current haplogroup to the ancestral root
  if(length(father) > 0) {
    path <- paste(father, find_path(father,relationships), sep = ">")
  }
  path <- gsub(">Root>Root", "", path)
  path <- unique(path)
  return(path)
}

# reverse path
find_and_reverse_path <- function(start_hap, relationships) {
  # use find_path function to find ancestral path
  path <- find_path(start_hap, relationships)
  final_path <- ifelse(grepl(start_hap, path), path, paste(start_hap, path, sep = ">"))
  
  # split,reverse,reassemble
  parts <- strsplit(final_path, split = ">", fixed = TRUE)[[1]]
  reversed_parts <- rev(parts)
  hap_path <- paste(reversed_parts, collapse = ">")
  
  return(list(final_path, hap_path))
}


### AGE estimate function
add_suffix <- function(time) {
  # Calculate the estimated age by subtracting the time from a reference year (1950) and rounding to the nearest decade
  result <- round(1950 - time, digits = -1)
  # If the result is negative, it's BCE; otherwise, it's CE. Format accordingly.
  if (result < 0) {
    return(paste(-result, "BCE", sep=" "))
  } else {
    return(paste(result, "CE", sep=" "))
  }
}

### Archaeological Era classification function
classify_archaeological_era <- function(mid_age) {
  # Check if the age is in BCE or CE and remove the suffix for further processing
  if (grepl(" BCE", mid_age)) {
    # BCE era
    mid_age <- gsub(" BCE", "", mid_age)
    mid_age <- as.numeric(mid_age)
    # Classify the age into historical periods based on the year
    if (mid_age >= 3300 && mid_age < 2000000) {
      return("Stone Age")
    } else if (mid_age >= 1200 && mid_age < 3300) {
      return("Bronze Age")
    } else if (mid_age >= 550 && mid_age < 1200) {
      return("Iron Age")
    } else {
      return("Unknown Era")
    }
  } else if (grepl(" CE", mid_age)) {
    # CE era
    mid_age <- gsub(" CE", "", mid_age)
    mid_age <- as.numeric(mid_age)
    
    if (mid_age >= 0 && mid_age < 476) {
      return("Imperial")
    } else if (mid_age >= 476 && mid_age < 1400) {
      return("Middle Ages")
    } else if (mid_age >= 1400 && mid_age < 1500) {
      return("Middle Ages/Modern")
    } else if (mid_age >= 1500) {
      return("Modern")
    } else {
      return("Unknown Era")
    }
  } else {
    return("Unknown Era")
  }
}

### Filter the dataframe to ensure the Full_Date column values are in increasing order
filter_increasing_full_date <- function(df) {
  df %>%
    filter(Full_Date >= cummax(lag(Full_Date, default = first(Full_Date) - 1)))
}

generate_final_table <- function(start_hap, country_name, data) {
  # Find hap path
  result <- find_and_reverse_path(start_hap, relationships)
  final_path <- result[[1]]
  rev_path <- result[[2]]
  
  # Transform the path into a dataframe for further processing
  path_table <- data.frame(Haplogroup = factor(unlist(strsplit(final_path, ">")), levels = unique(unlist(strsplit(final_path, ">")))), stringsAsFactors = FALSE)
  
  # Merge the path table with the data, filtering by country if specified
  if (is.null(country_name) || country_name == "" || country_name == "Default") {
    matched_data <- merge(path_table, data, by.x="Haplogroup", by.y = "Y_haplo", all.x = TRUE)
  } else {
    filtered_data <- data[data$Political_Entity == country_name,]
    matched_data <- merge(path_table, filtered_data, by.x="Haplogroup", by.y = "Y_haplo", all.x = TRUE)
  }
  
  # Delete NA rows
  new_data <- na.omit(matched_data)
  new_data <- new_data %>%
    group_by(Haplogroup) %>%
    arrange(Full_Date, .by_group = TRUE) %>%
    ungroup()
  new_data <- new_data %>% distinct(Haplogroup, .keep_all = TRUE)
  
  new_data <- filter_increasing_full_date(new_data)
  new_data$Age_Estimate <- sapply(new_data$Full_Date, add_suffix)
  new_data$Archaeological_Era <- sapply(new_data$Age_Estimate, classify_archaeological_era)
  
  ### Calculate time passed
  # diff: calculates the difference between adjacent elements in a vector
  new_data <- new_data %>%
    mutate(Time_Passed = c(round(diff(Full_Date), -1), NA)) %>%
    mutate(Time_Passed = if_else(row_number() == n(), "", paste0(Time_Passed, " years")))
  
  # Finalize
  final_table <- subset(new_data, select = -c(Political_Entity, Full_Date, country_frequency, type_haps,
                                              total_num_haps, country_population, global_frequency,
                                              global_population))
  
  if (country_name == "Default") {
    final_table <- subset(final_table, select = -c(country_count, country_individuals))
    final_table <- select(final_table, Haplogroup, Age_Estimate, Archaeological_Era,
                          Time_Passed, global_count, global_individuals)
    names(final_table) <- c("Haplogroup", "Age Estimate", "Archaeological Era", 
                            "Time Passed", "Immediate Desendants", "Tested Modern Descendants")
  } else { 
    final_table <- subset(final_table, select = -c(global_count, global_individuals))
    final_table <- select(final_table, Haplogroup, Age_Estimate, Archaeological_Era,
                          Time_Passed, country_count, country_individuals)
    names(final_table) <- c("Haplogroup", "Age Estimate", "Archaeological Era", 
                            "Time Passed", "Immediate Desendants", "Tested Modern Descendants")
  }
  return(list(path = rev_path, table = final_table))
}

# generate_final_table(merged_data)



