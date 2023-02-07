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
from tokenizers import Tokenizer
from tokenizers.models import BPE
import io, csv
## Some trial
fastagz_file = 'data/version2.5.2019-10-09/An-1.chr.all.v2.0.fasta.gz'


def parse_raw_data(raw_data_file):
    '''
    Read raw .fasta file.

    Arguments:
    raw_data_path -- str, path to data file.

    Returns:
    chromosome_list -- All the sequence put together
    '''
    chromosome_list = []

    with gzip.open(raw_data_file, "rt") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            chromosome_list.append(str(record.seq).replace('N', '').upper()) # Just omit N as if it were not there and put all uper case
            if 'chr' not in record.id:
                continue
    return ''.join(chromosome_list)

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

def output_data(kmer_sequence: str):
    with open('./data/train.txt', 'w') as f:
        f.write(kmer_sequence)


def BPE():
    pass

my_sequence = parse_raw_data(fastagz_file)


output_data(seq2kmer(my_sequence, 5))