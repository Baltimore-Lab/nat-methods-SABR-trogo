#!/usr/bin/env python
# 
# Created by: Michael Leonard
# Last updated: November 26, 2018
#
# Usage:   mutant_analysis.py sample.sam epitope_start epitope_stop library.txt
# Example: mutant_analysis.py sample.sam 5292 5327 library.txt
# Output:  stdout 
#
# Extracts nucleotide sequences encoding specific epitopes from a .sam alignment
# to a single-chain-trimer, converts to amino acid, and then counts the number 
# of hits from a list of provided epitopes in a library
#
# "epitope_start" and "epitope_stop" refer to genomic position of the epitope 
#  in the aligned plasmid. Script is not CIGAR string aware.
#
# 
# Prerequisite analysis:
#   Align data:
#      REF="/Reference/A2-SABR/"
#      bwa mem -B 1 "$REF"/A2 sample.fastq > sample.sam
#        
#   Analyze data:
#      mutant_analysis.py sample.sam 5292 5327 library.txt > sample.counts.txt
#    
#

import sys

filename = sys.argv[1]
start_position = int(sys.argv[2])
stop_position = int(sys.argv[3])
library_filename = sys.argv[4]
target_length = stop_position + 1 - start_position

mutant_count_table = {}
mutant_count_table["Unclassified"] = 0

codon_table = {
'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
'TAC':'Y', 'TAT':'Y', 'TAA':'.', 'TAG':'.',
'TGC':'C', 'TGT':'C', 'TGA':'.', 'TGG':'W',
}

sys.stderr.write("Reading library file\n")
sys.stderr.flush()

# Generate an ordered list of the library
# And fill dictionary at the same time
mutant_library_file = open(library_filename, 'r')
mutant_library = []
for mutant in mutant_library_file:
    mutant = mutant.strip()
    mutant_library.append(mutant)
    mutant_count_table[mutant] = 0
mutant_library_file.close()


# Process .sam file
sys.stderr.write("Reading sam file\n")
sys.stderr.flush()

samfile = open(filename, 'r')
read_count = 0
line = samfile.readline()

while line:
    # Show progress every 100,000 lines
    read_count += 1
    if read_count % 100000 == 0:
        sys.stderr.write("Processed %s reads\n" % read_count)
        sys.stderr.flush()
    
    # Skip comment lines
    if line[0] == "@":
        line = samfile.readline()
        continue
    
    # Parse sam format
    line_list = line.split('\t')
    read_name = line_list[0]
    read_flag = int(line_list[1])
    read_reference = line_list[2]
    read_start = int(line_list[3])
    read_mapping_quality = int(line_list[4])
    read_cigar = line_list[5]
    read_pair_name = line_list[6]
    read_pair_start = int(line_list[7])
    read_insert_length = int(line_list[8])
    read_sequence = line_list[9]
    read_quality = line_list[10]
    read_comments = line_list[11:]
    
    read_length = len(read_sequence)
    
    # Skip any alignment that is not forward direction
    # and also skip secondary alignments
    if read_flag != 0:
        line = samfile.readline()
        continue
    
    # Extract target sequence
    # Not CIGAR aware
    if not (start_position >= read_start and stop_position <= (read_start + read_length)):
        line = samfile.readline()
        continue
        
    target_sequence_start = start_position - read_start
    target_sequence_stop = stop_position - read_start + 1
    target_sequence = read_sequence[target_sequence_start:target_sequence_stop]
    
    # Skip target sequences with N
    if 'N' in target_sequence:
        line = samfile.readline()
        continue
    
    # DNA to AA translation 
    amino_acid_sequence = ""    
    for amino_acid in range(target_length/3):
        amino_acid_start = amino_acid*3
        amino_acid_stop = amino_acid*3 + 3
        codon = target_sequence[amino_acid_start:amino_acid_stop]
        if len(codon) == 3:
            amino_acid_sequence += codon_table[codon]
        else:
            sys.stderr.write("AA seq truncated: Read %s, Target seq: %s\n" % (read_name, target_sequence))
            line = samfile.readline()
            continue
    
    # Add epitope to table
    if amino_acid_sequence in mutant_count_table:
        mutant_count_table[amino_acid_sequence] += 1
    else:
        mutant_count_table["Unclassified"] += 1
    
    line = samfile.readline()
    
samfile.close()

#for key in mutant_count_table:
#    print("%s\t%s" % (key, mutant_count_table[key]))


# Output counts in order of library file
sys.stderr.write("Outputting final counts\n")
sys.stderr.flush()

for mutant in mutant_library:
    if mutant in mutant_count_table.keys():
        print("%s\t%s" % (mutant, mutant_count_table[mutant]))
    else:
        print("%s\t%s" % (mutant, 0))
        

print("%s\t%s" % ("Unclassified", mutant_count_table["Unclassified"]))

