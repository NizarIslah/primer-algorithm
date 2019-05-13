"""
Main file, run the program from here.
"""

import time
from PrimerDesign import *

if __name__ == '__main__':
    genome_fa = "GCF_000250985.1_Nema_parisii_ERTm1_V3_rna.fna"
    prot_table = "ProteinTableTrim.txt"
    genome_size = 2660  # choose how many of the genes in the genome you want to run the algorithm on
    size_min = 25  # minimum pool size
    size_max = 100  # max pool size
    step = 25
    tm_cutoff = 55
    kmer_file = "DNA_20mers.txt"  # file containing randomly generated 20 nucleotide DNA sequences needed for the oligos
    out_file = "results"
    fmt = ".csv"

    # Start the timer
    t0 = time.time()

    # run the algorithm on different pool sizes
    for size in range(size_min, size_max + 1, step):
        # extract gene information from the fasta file and protein table (can specify file names as parameters)
        primers, genes, pos = load_genes(genome_size, genome_fa, prot_table)

        # Run the pooling algorithm to place primers in pools based on specific restraints including size, dist, TM.
        pool = pool_placement(size, tm_cutoff, primers, genes, pos)

        # Create the final oligos for each set of primers, required for initial round of PCR.
        oligos = create_oligos(pool, kmer_file)

        # save the sequences to an output file
        save_output(out_file, fmt, size, oligos)

    print("Runtime (seconds):", time.time() - t0)
