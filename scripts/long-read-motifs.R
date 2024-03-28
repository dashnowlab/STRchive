lr100_allele_lengths = read.csv('/Users/quinlan/Downloads/all_samples.allele_lengths.tsv', sep = '\t', stringsAsFactors = FALSE)
lr100_motif_counts = read.csv('/Users/quinlan/Downloads/all_samples.motif_counts.tsv', sep = '\t', stringsAsFactors = FALSE)

library(dplyr)
library(tidyr)

# Convert MotifCount column to numeric
lr100_motif_counts$MotifCount <- as.numeric(lr100_motif_counts$MotifCount)

lr100_motif_counts <- lr100_motif_counts %>%
  filter(MotifCount != 0)

# Group lr100_motif_counts by Locus and Motif, then count the occurrences of each unique Motif
unique_motif_counts <- lr100_motif_counts %>%
  group_by(Locus, Motif) %>%
  summarise(TotalMotifCount = sum(MotifCount)) %>%
  ungroup() %>%
  group_by(Locus) %>%
  mutate(UniqueMotifCount = n()) %>%
  select(Locus, Motif, TotalMotifCount)  # Optional: Select only necessary columns

# Create a dataframe with unique motifs as columns
motif_counts_pivoted <- unique_motif_counts %>%
  pivot_wider(names_from = Motif, values_from = TotalMotifCount, values_fill = 0)

# Sum the number of non-zero values for each row, excluding the first column
motif_counts_pivoted$UniqueMotifCount <- rowSums(motif_counts_pivoted[, -1] > 0)

motif_unique_counts <- motif_counts_pivoted %>%
  select(Locus, UniqueMotifCount) %>%
  filter(UniqueMotifCount != 1)

### DMD analysis
lr100_DMD_path_counts <- lr100_motif_counts %>%
  filter(Locus == "DMD_DMD") %>%
  filter(MotifCount >= 59)
