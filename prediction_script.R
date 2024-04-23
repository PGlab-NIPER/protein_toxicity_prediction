#Call and install the dependencies 
library("caret")
library("protr")
library("Boruta")
library("ranger")
library("readr")

# Input: User can upload the query sequences in the form of a FASTA file that can contain n number of sequences.
# Make sure to remove the ambiguous amino acids (such as X) from the sequences; otherwise, it will result in an error.
fasta_file <- "testfile.fasta"

# Alternatively, users can choose the FASTA file from their directory using the file.choose() function.
# This allows for more flexibility in selecting the input file.
#fasta_file <- file.choose()

# Each sequence is stored as a list element in the 'sequences' variable
sequences <- readFASTA(fasta_file)

# Create a list to store the amino acid composition based features for each sequence
composition_feature_list <-list()
# Iterate over each sequence in the 'sequences' variable
for (i in seq_along(sequences)) {
  # Load the current sequence
  X <- sequences[[i]]
  # Calculate the amino acid composition based features
  aa_composition <- extractAAC(X)
  dpc_composition <- extractDC(X)
  tpc_composition <- extractTC(X)
  ctriad_feature <- extractCTriad(X)
  compo <- extractCTDC(X)
  distri <- extractCTDD(X)
  transi <- extractCTDT(X)
  # Convert the features to data frames
  aac <- as.data.frame(aa_composition)
  dpc <- as.data.frame(dpc_composition)
  tpc <- as.data.frame(tpc_composition)
  ctriad <- as.data.frame(ctriad_feature)
  composition <- as.data.frame(compo)
  transition <- as.data.frame(transi)
  distribution <- as.data.frame(distri)
  # Rename the columns of each data frame
  colnames(aac) <- 'Column'
  colnames(dpc) <- 'Column'
  colnames(tpc) <- 'Column'
  colnames(ctriad) <- 'Column'
  colnames(composition) <- 'Column'
  colnames(transition) <- 'Column'
  colnames(distribution) <- 'Column'
  # Combine all features into a single data frame
  XX <- rbind(aac, dpc, tpc, ctriad, composition, transition, distribution)
  XX <- t(XX)
  # Add the features to the composition feature list
  composition_feature_list[[i]] <- XX   
}

# Convert list of dataframes to a single dataframe
composition_df <- do.call(rbind, composition_feature_list)
# Set row names to NULL
rownames(composition_df) <- NULL

# Initialize an empty list to store PnGT features
PnGT_feature_list <- list()

# The PnGT (Protein n-gram Tool) calculates features based on predefined amino acid patterns.
# These patterns represent specific combinations of amino acids and their properties.
# To mimic the PnGT tool's functionality, we read a file named "Pattern.xlsx", which contains
# these predefined patterns along with their corresponding feature values.
# By using these patterns, we can efficiently calculate PnGT-based features for protein sequences
# by matching the amino acid combinations in the sequences with the predefined patterns.
Pattern <- read_excel("Pattern.xlsx")

# Loop through each sequence in the list of sequences
for (i in seq_along(sequences)) {
  #get the current sequence
  X <- sequences[[i]]
  # Calculate PnGT based features
  text <- X
  chars <- strsplit(text, "")[[1]]
  # Generate all possible two-letter combinations of amino acids
  combinations <- vector("list", length = length(chars) - 1)
  for (j in 2:length(chars)) {
    combinations[[j - 1]] <- paste(chars[j - 1], chars[j], sep = "")
  }
  combinations <- unlist(combinations)
  # Initialize a dataframe to store PnGT features for the current sequence
  aa <- matrix(NA, ncol = 12, nrow = 1)
  aa <- as.data.frame(aa)
  colnames(aa) <- c('Seq', c(1:11))
  aa[, 2:12] <- 0
  # Iterate over each two-letter combination of amino acids
  for (j in 1:length(combinations)) {
    text_to_check <- combinations[j]
    # Compare the current combination with the predefined patterns in Pattern.xlsx
    for (k in 1:nrow(Pattern)) {
      possible_characters <- as.character(Pattern[k, 2])
      is_present <- unlist(strsplit(text_to_check, '')) %in% unlist(strsplit(possible_characters, ','))
      # If the current combination matches a predefined pattern, increment the corresponding feature count
      sum1 <- sum(as.numeric(is_present))
      if (sum1 == 2) {
        aa[1, k + 1] <- aa[1, k + 1] + 1
      }
    }
  }
  # Add the sequence to the dataframe
  aa$Seq <- text
  aa <- aa[, -1]
  colnames(aa) <- c('Tiny', 'Small', 'Aliphatic', 'Nonpolar', 'Aromatic', 'Polar', 'Charged', 'Basic', 'Acidic', 'Hydrophobic', 'Hydrophilic')
  # Store the PnGT features for the current sequence in the PnGT_feature_list
  PnGT_feature_list[[i]] <- aa
}

# Combine the list of PnGT features into a single datafram
PnGT_feature_df <- do.call(rbind, PnGT_feature_list)
# Set row names to NULL
rownames(PnGT_feature_df) <- NULL

##Now combine the both type of feature together into one file
# Create empty lists to store data frames
composition_dfs <- list()
PnGT_dfs <- list()
# Extract data frames from composition_feature_list and PnGT_feature_list
for (i in seq_along(composition_feature_list)) {
  composition_dfs[[i]] <- composition_feature_list[[i]]
  PnGT_dfs[[i]] <- PnGT_feature_list[[i]]
}
# Combine data frames row-wise
composition_df <- do.call(rbind, composition_dfs)
PnGT_df <- do.call(rbind, PnGT_dfs)
# Combine composition_df and PnGT_df column-wise
combined_df <- cbind(composition_df, PnGT_df)

#Load all three best performing models saved as pickle in pickle folder
#Load model saved as pickle
RF1 <- readRDS("pickle/best_perf_cardio_model_RF.pickle")
RF2 <- readRDS("pickle/best_perf_neuro_model_RF.pickle")
SVM1 <- readRDS("pickle/best_perf_entero_model_SVM.pickle")

# Initialize an empty data frame to store results
results <- list()
# Iterate over each row of the combined data frame
for (i in 1:nrow(combined_df)) {
  # Get the features for the current sequence
  features <- combined_df[i, ]
  which(colnames(features)=="NA")
  colnames(features)[23]="NA."
  # Predict toxicity using the three models
  toxpred1 <- predict(RF1, newdata = features)
  toxpred2 <- predict(RF2, newdata = features)
  toxpred3 <- predict(SVM1, newdata = features)
  
  # Determine toxicity terms based on predictions
  term1 <- ifelse(toxpred1 == 1, 'Cardiotoxic', 'Non-Cardiotoxic')
  term2 <- ifelse(toxpred2 == 1, 'Neurotoxic', 'Non-Neurotoxic')
  term3 <- ifelse(toxpred3 == 1, 'Enterotoxic', 'Non-Enterotoxic')
  
  # Create a data frame for the current sequence and its predictions
  result <- data.frame(Sequence = sequences[[i]], Cardio = term1, Neuro = term2, Entero = term3)
  
  # Append the result to the results list
  results[[i]] <- result
}
# Combine all results into a single data frame
final_results <- do.call(rbind, results)
# Write results to CSV file
write.csv(final_results, file = 'Result.csv', row.names = FALSE)


