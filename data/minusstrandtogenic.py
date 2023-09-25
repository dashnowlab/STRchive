import pandas as pd
from Bio.Seq import Seq

# Read the CSV file into a DataFrame
df = pd.read_csv("STR-disease-loci.csv")

# Check if the gene_strand is negative and the new pathogenic_motif_reference_orientation
# matches repeatunit_pathogenic_geneorientation before updating
negative_strand = (df['gene_strand'] == '-') & (df['pathogenic_motif_reference_orientation'] == df['repeatunit_pathogenic_geneorientation'])

# Update the pathogenic_motif_gene_orientation column for matching rows
df.loc[negative_strand, 'pathogenic_motif_gene_orientation'] = df.loc[negative_strand, 'pathogenic_motif_reference_orientation'].apply(lambda x: str(Seq(x).reverse_complement()))

# Save the modified DataFrame back to the CSV file
df.to_csv("STR-disease-loci_new.csv", index=False)

print("Done! Check STR-disease-loci.csv for the updated values.")