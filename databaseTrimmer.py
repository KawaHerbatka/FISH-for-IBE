#This script can be used to filter out sequences from a fasta file, based on another file, that contains
#just the names of sequences to filter out. Output is a new fasta file without specified sequences
#This script was made for use with oligoNdesign as part of filtering out target sequences from database.

import argparse
import re

###########################################################################################################
#Input arguments 

parser = argparse.ArgumentParser(prog="BlastSeeker", description="Placeholder", epilog="")

parser.add_argument('-o', action="store", required=True, dest='out_folder', help="place to save output as .fasta file")

parser.add_argument('-x', action="store", required=True, dest='cut_out', help="file with sequences to be removed")

parser.add_argument('-d', action="store", required=True, dest='db_file', help="database file as .fasta file")

args = parser.parse_args()

#############################################################################################################

with open(args.cut_out, 'r') as cut_file, open(args.db_file, 'r') as db_file, open(args.out_folder+"/modified_db.fasta", 'w') as out_folder:
    cutNext = False
    for line in db_file:
        cut_file.seek(0)
        for seq_id in cut_file:
            if re.sub('\n','',seq_id) in line:
                cutNext = True
        if not cutNext:
            print(line, end="", file=out_folder)
        elif '>' not in line:
            cutNext = False

        








