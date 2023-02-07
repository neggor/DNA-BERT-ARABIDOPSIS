# Ok so I have 7 whole sequence genomes from arabidopsis, and I want to use this as my pre-training data!
from Bio import SeqIO
from pathlib import Path

# Automate data downloading
import sys
import os
import wget
import gzip
import shutil
import subprocess

 
## Some trial
fastq_folder = 'data/version2.5.2019-10-09/An-1.chr.all.v2.0.fasta'


def parse_raw_data(raw_data_path):
    '''
    Read raw .fasta file.

    Arguments:
    raw_data_path -- str, path to data file.

    Returns:
    chromosome_list -- list of str, one entry per chromosome.
    '''
    chromosome_list = []
    for record in SeqIO.parse(raw_data_path, "fasta"):
        chromosome_list.append(str(record.seq))
    
    return chromosome_list

def seq2kmer(seq, k):
    """
    Convert original sequence to kmers
    
    Arguments:
    seq -- str, original sequence.
    k -- int, kmer of length k specified.
    
    Returns:
    kmers -- str, kmers separated by space
    """
    kmer = [seq[x:x+k] for x in range(len(seq)+1-k)]
    kmers = " ".join(kmer)
    return kmers

my_chromosome_list = parse_raw_data(fastq_folder)

#print(seq2kmer(my_chromosome_list[0], 5))