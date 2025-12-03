install.packages("baseq")
library(baseq)

file_path = "/media/sf_BioMachine_files/Euglenibacteraceae/anti-eugb/6th probe_mismatches.fasta"
out = "/media/sf_BioMachine_files/Euglenibacteraceae/anti-eugb"
# Write to custom directory
write.rna_to_dna(file_path, output_dir = out)
