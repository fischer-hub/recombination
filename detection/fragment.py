#!/usr/bin/env python3

import argparse
import subprocess
import csv

# read in arguments
parser = argparse.ArgumentParser(description="Takes a FASTA file as input and returns a multi-FASTA file with every sequence being a fragment of the input sequence of size WINDOW_SIZE with offset STEP_SIZE to the previous fragment.", epilog="Example usage: python3 fragment.py <FASTA_FILE> <WINDOW_SIZE> <STEP_SIZE>" )
parser.add_argument("FASTA_FILE", help="FASTA file containing the input sequence.")
parser.add_argument("WINDOW_SIZE", help="Size of each fragment.")
parser.add_argument("STEP_SIZE", help="Size of shift of the window between fragments.")
parser.add_argument("--full", help="Also run pangolin on created fragments.")
args = parser.parse_args()

window_sz    = int(args.WINDOW_SIZE)
step_sz      = int(args.STEP_SIZE)
counter      = 0
in_filename  = args.FASTA_FILE
out_filename = f"{in_filename.rsplit('.', 1)[0]}_fragments_{window_sz}_{step_sz}.fasta"


# read in input sequence
with open(args.FASTA_FILE, 'r') as fasta:
    #next(fasta)
    seq_id = fasta.readline()
    sequence = fasta.read()


# write fragment multi-FASTA file
with open(f"fragLEARN/{out_filename}", 'w') as fragments:
    for i in range(0, len(sequence)-window_sz+1, step_sz):
        counter += 1
        print(f"{seq_id.rstrip()}_fragment:{counter}_startpos:{i}_endpos:{i+window_sz-1}\n{sequence[i:i+window_sz]}", file = fragments)
        #print(f">fragment:{counter}_startpos:{i}_endpos:{i+window_sz-1}\n{'N'*i}{sequence[i:i+window_sz]}{'N'*(len(sequence)-i+window_sz)}", file = fragments) #masked

counter = 0
with open(f"fragLEARN/metadata.csv", 'a') as meta:
    print("sequence_name,lineage", file = meta)
    for i in range(0, len(sequence)-window_sz+1, step_sz):
        counter += 1
        print(f"{seq_id[1:].rstrip()}_fragment:{counter}_startpos:{i}_endpos:{i+window_sz-1},{in_filename.rsplit('.', 1)[0]}", file = meta)

if args.full:
    # run pangolin on multi-FASTA file
    subprocess.run(["pangolin", out_filename, "--outfile", f"{out_filename}_pangolin.csv", "--min-length", f"{window_sz}"], check = True) #, "--max-ambig", "0.99"]


    # check pangolin output for classifications
    with open(f"{out_filename}_pangolin.csv", 'r') as lin_report:
        reader = csv.reader(lin_report, delimiter=",", quotechar='"')
        next(reader, None)
        data_read = [row for row in reader]

    dd = dict()

    # evaluate lineage classifications from pangolin output
    for row in data_read:
        dd[row[1]] = 1 if row[1] not in dd else dd[row[1]]+1

    with open(f"fragment_classification_{window_sz}_{step_sz}.txt", 'a') as out:
        print(f">{in_filename.rsplit('.', 1)[0]}", file = out)
        for lineage, count in dd.items() :
            print(f"lineage: {lineage: <12} classsifications {count: <4} ({(count/len(data_read)) * 100} %)", file = out)
