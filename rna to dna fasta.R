install.packages("baseq")
library(baseq)

file_path = "/path/to/rna/fast/file.fasta"
out = "/path/to/output/FOLDER"
# Write to custom directory
write.rna_to_dna(file_path, output_dir = out)
