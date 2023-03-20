# -*- coding: utf-8 -*-
"""
Created on Fri Mar 17 10:35:30 2023

@author: jacki
"""

from Bio import Entrez
from Bio import SeqIO
import os
import numpy as np

#--STEP 2--#
Entrez.email = ''
# Retrieve the GenBank record
handle = Entrez.efetch(db='nucleotide', id='NC_006273.2', rettype='gb', retmode='text')
hcmv_record = handle.read()
handle.close()
print(hcmv_record)
hcmv_gb = SeqIO.read(hcmv_record, 'genbank')
# Extract the CDS features from the record
cds_features = [f for f in hcmv_gb.features if f.type == 'CDS']
# Create a dictionary with RefSeq protein_id as key and CDS sequence as value
cds_seqs = {}
for cds in cds_features:
    # Extract the RefSeq protein_id from the CDS feature
    protein_id = [qualifier.value for qualifier in cds.qualifiers['protein_id'] if 'RefSeq' in qualifier.description][0]
    # Extract the CDS sequence from the feature and remove any gaps or stop codons
    cds_seq = str(cds.extract(hcmv_gb.seq)).replace('-', '').rstrip('N')
    # Store the CDS sequence in the dictionary
    cds_seqs[protein_id] = cds_seq
# Write the CDS sequences to a FASTA file
with open('hcmv_cds.fasta', 'w') as f:
    for protein_id, cds_seq in cds_seqs.items():
        f.write(f'>{protein_id}\n{cds_seq}\n')
# Log the number of CDS features
num_cds = len(cds_features)
print(f'The HCMV genome (NC_006273.2) has {num_cds} CDS.')

#--STEP 3--#
# Define the paths to the transcriptome index and the FASTQ files for each sample
transcriptome_index = 'testindex.idx'
sample_files = {'sample1': ('SRX2896360_1.fastq', 'SRX2896360_2.fastq'),
                'sample2': ('SRX2896363_1.fastq', 'SRX2896363_2.fastq'),
                'sample3': ('SRX2896374_1.fastq', 'SRX2896374_2.fastq'),
                'sample4': ('SRX2896375_1.fastq', 'SRX2896375_2.fastq')}
# Loop through each sample and run kallisto
for sample, (fastq1, fastq2) in sample_files.items():
    output_dir = f'{sample}_output'
    os.makedirs(output_dir, exist_ok=True)
    cmd = f'kallisto quant -i {transcriptome_index} -o {output_dir} -b 10 {fastq1} {fastq2}'
    os.system(cmd)  
# Define the paths to each sample's output directory
sample_dirs = ['sample1_output', 'sample2_output', 'sample3_output', 'sample4_output']
# Calculate TPM statistics for a single sample
def calculate_tpm_stats(sample_dir):
    # Read in the abundance.tsv file
    with open(os.path.join(sample_dir, 'abundance.tsv'), 'r') as f:
        lines = f.readlines()
    # Extract the CDS ids from the header row of the kallisto output file
    cds_ids = lines[0].strip().split('\t')[4:]
    # Parse the lines to create a dictionary of CDS TPM values
    cds_tpm = {line.split('\t')[0]: float(line.split('\t')[4]) for line in lines}
    # Extract the TPM values for each CDS from the kallisto output file
    tpm_values = [list(map(float, line.strip().split('\t')[4:])) for line in lines[1:]]
    # Calculate the minimum, median, mean, and maximum TPM values for the CDSs
    min_tpm = [min(tpm_values[i]) for i in range(len(cds_ids))]
    med_tpm = [np.median(tpm_values[i]) for i in range(len(cds_ids))]
    mean_tpm = [np.mean(tpm_values[i]) for i in range(len(cds_ids))]
    max_tpm = [max(tpm_values[i]) for i in range(len(cds_ids))]
    
#--STEP 5--#    
# Creation of Blast DB   
Entrez.email = ''
# Search for nucleotide sequences in the Betaherpesvirinae subfamily
handle = Entrez.esearch(db="nucleotide", term="Betaherpesvirinae[subtree]", retmax=100000)
record = Entrez.read(handle)
# Download in Fasta format
handle = Entrez.efetch(db="nucleotide", id=record["IdList"], rettype="fasta", retmode="text")
sequences = handle.read()
# Write the sequences to a FASTA file
with open("betaherpesvirinae.fasta", "w") as outfile:
    outfile.write(sequences)
# Create a blast database using the FASTA file
os.system("makeblastdb -in betaherpesvirinae.fasta -dbtype nucl")
# Set up the BLAST search
query_seq = SeqIO.read("most_diff_expr_prot.fasta", "fasta")
result_handle = NCBIWWW.qblast("blastx", "nt", query_seq.seq, entrez_query="Betaherpesvirinae[subtree]", hitlist_size=10)
# Parse the BLAST results
blast_records = NCBIXML.parse(result_handle)
for record in blast_records:
    # Write the header row
    with open("blast_results.txt", "w") as outfile:
        outfile.write("sacc\tpident\tlength\tqstart\tqend\tsstart\tsend\tbitscore\tevalue\tstitle\n")
    # Write the hits to the log file
    for alignment in record.alignments:
        hit = alignment.hsps[0]
        with open("blast_results.txt", "a") as outfile:
            outfile.write(f"{alignment.accession}\t{hit.identities/alignment.length*100:.2f}\t{hit.align_length}\t{hit.query_start}\t{hit.query_end}\t{hit.sbjct_start}\t{hit.sbjct_end}\t{hit.bits}\t{hit.expect}\t{alignment.title}\n")
# Read the blast results into a pandas dataframe
blast_results = pd.read_csv("blast_results.txt", sep="\t")
# Group the hits by subject accession and summarize the top hits
top_hits = blast_results.groupby("sacc").agg({
    "pident": "mean",
    "length": "mean",
    "bitscore": "max",
    "evalue": "min",
    "stitle": "first"
}).reset_index()
# Write the summary to the log file
with open("blast_summary.txt", "w") as outfile:
    outfile.write("sacc\tpident_mean\tlength_mean\tbitscore_max\tevalue_min\tstitle_first\n")
    for _, row in top_hits.iterrows():
        outfile.write(f"{row['sacc']}\t{row['pident']:.2f}\t{row['length']:.0f}\t{row['bitscore']:.0f}\t{row['evalue']:.2e}\t{row['stitle']}\n")

