#!/usr/bin/env python3

#This script takes blast results in, cuts out all under given percentage of identical positions and 
#saves a .csv file with sequence id in DB, perc. of identity and sequence itself.
#It was used for FISH porbe design with oligoNdesign to filter out target sequence from SILVA.

import re
import csv
import argparse

###########################################################################################################
#Input arguments 

parser = argparse.ArgumentParser(prog="BlastSeeker", description="Placeholder", epilog="")

parser.add_argument('-o', action="store", required=True, dest='out_folder', help="place to save output as .csv file")

parser.add_argument('-b', action="store", required=True, dest='blast_out', help="file with blast results")

parser.add_argument('-d', action="store", required=True, dest='db_file', help="database file as .fasta file")

args = parser.parse_args()

############################################################################################################
#Read the blast output data

data = {}

with open(args.blast_out, 'r') as file:
    for line in file:
        if float(line.split("\t")[2]) >= 90:
            if line.split("\t")[1] not in data:
                data.update({line.split("\t")[1]:[line.split("\t")[2], ""]})

###########################################################################################################
#Get the sequence itself from the database

with open(args.db_file, 'r') as file:
    for id in data:
        file.seek(0)
        takeNext = False
        for line in file:
            if id in line:
                takeNext = True
            elif takeNext:
                data[id][1] = re.sub("\n", "", line)
                break

##########################################################################################################
#save the .csv file with all the results

with open(args.out_folder+"/csv_file.csv", 'w') as csvfile:
    csvwriter = csv.writer(csvfile)
    csvwriter.writerow(['Sequence id', 'Perc. Identity', 'Sequence'])
    for id in data:
        csvwriter.writerow([id, data[id][0], data[id][1]])
        


