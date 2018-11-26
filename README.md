# nat-methods-SABR-trogo
Code for analysis of SABR and Trogocytosis antigen discovery assay sequencing

## mutant_analysis.py
Created by: Michael Leonard  
Last updated: November 26, 2018  

Usage:   mutant_analysis.py sample.sam epitope_start epitope_stop library.txt  
Example: mutant_analysis.py sample.sam 5292 5327 library.txt  
Output:  stdout 

Extracts nucleotide sequences encoding specific epitopes from a .sam alignment
to a single-chain-trimer, converts to amino acid, and then counts the number
of hits from a list of provided epitopes in a library.

"epitope_start" and "epitope_stop" refer to genomic position of the epitope 
in the aligned plasmid. Script is not CIGAR string aware.

  
Prerequisite analysis:  
Align data:  
REF="/Reference/A2-SABR/"  
bwa mem -B 1 "$REF"/A2 sample.fastq > sample.sam  
