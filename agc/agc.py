#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""OTU clustering"""

import argparse
import sys
import os
import gzip
import statistics
import textwrap
from collections import Counter
# https://github.com/briney/nwalign3
# ftp://ftp.ncbi.nih.gov/blast/matrices/
import nwalign3 as nw

__author__ = "Jeyani George Clement"
__copyright__ = "Universite Paris Cité"
__credits__ = ["Jeyani George Clement"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Jeyani George Clement"
__email__ = "jeyanigeorgeclement@gmail.com"
__status__ = "Developpement"


def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', '-amplicon_file', dest='amplicon_file', type=isfile, required=True, 
                        help="Amplicon is a compressed fasta file (.fasta.gz)")
    parser.add_argument('-s', '-minseqlen', dest='minseqlen', type=int, default = 400,
                        help="Minimum sequence length for dereplication (default 400)")
    parser.add_argument('-m', '-mincount', dest='mincount', type=int, default = 10,
                        help="Minimum count for dereplication  (default 10)")
    parser.add_argument('-c', '-chunk_size', dest='chunk_size', type=int, default = 100,
                        help="Chunk size for dereplication  (default 100)")
    parser.add_argument('-k', '-kmer_size', dest='kmer_size', type=int, default = 8,
                        help="kmer size for dereplication  (default 10)")
    parser.add_argument('-o', '-output_file', dest='output_file', type=str,
                        default="OTU.fasta", help="Output file")
    return parser.parse_args()

def read_fasta(amplicon_file, minseqlen):
    seq = ""
    with gzip.open(amplicon_file, "rt") as filin:
        for line in filin:
            if line.startswith(">"):
                # dans le cas où on a déjà stocké une seq et qu'on veut pas la perdre
                if len(seq) >= minseqlen:
                    yield seq
                seq = ""
            else:
                seq += line.strip()
        if len(seq) >= minseqlen:
            yield seq

def dereplication_fulllength(amplicon_file, minseqlen, mincount):
    seq = list(read_fasta(amplicon_file, minseqlen))
    seq_unique = set(seq)
    count = 0
    seq_count = []

    for sequence in seq_unique:
        if sequence in seq:
            count = seq.count(sequence)
            if count >= mincount:
                seq_count.append([sequence, count])
    
    seq_count.sort(key = lambda x:x[1], reverse = True)
    
    for seq in seq_count:
        yield seq


def get_identity(alignment_list):
    """Prend en une liste de séquences alignées au format ["SE-QUENCE1", "SE-QUENCE2"]
    Retourne le pourcentage d'identite entre les deux."""
    seq_1 = alignment_list[0]
    seq_2 = alignment_list[1]

    common = 0
    # calcul du nombre de nucléotides identiques 
    for i in range(len(seq_1)):
        if seq_1[i] == seq_2[i]:
            common += 1
    
    id = (common / len(seq_1))*100
    return id

def abundance_greedy_clustering(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    seq = list(dereplication_fulllength(amplicon_file, minseqlen, mincount))
    otu_list = []
    otu_list.append(seq[0])

    for sequence in seq[1:]:
        for otu in otu_list:
            align_list = nw.global_align(sequence[0], otu[0], gap_open=-1, gap_extend=-1, 
            matrix=os.path.abspath(os.path.join(os.path.dirname(__file__), "MATCH")))
            id = get_identity(align_list)

            if id < 97:
                otu_list.append(sequence)

        return otu_list
            
def write_OTU(OTU_list, output_file):
    with open(output_file, "w") as filout: 
        for i in range(len(OTU_list)):
            filout.write(f">OTU_{i+1} occurrence:{OTU_list[i][1]}\n")
            filout.write(f"{textwrap.fill(OTU_list[i][0], width=80)}\n")


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    # Votre programme ici
    sequence = read_fasta(args.amplicon_file, args.minseqlen)

#==============================================================
# Chimera removal section
#==============================================================

def get_unique(ids):
    return {}.fromkeys(ids).keys()

def common(lst1, lst2): 
    return list(set(lst1) & set(lst2))

def get_chunks(sequence, chunk_size):
    """Split sequences in a least 4 chunks
    """
    pass

def cut_kmer(sequence, kmer_size):
    """Cut sequence into kmers"""
    pass

def get_unique_kmer(kmer_dict, sequence, id_seq, kmer_size):
    pass

def detect_chimera(perc_identity_matrix):
    pass

def search_mates(kmer_dict, sequence, kmer_size):
    pass

def chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    pass


if __name__ == '__main__':
    main()
