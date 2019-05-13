# primer-algorithm

## Overview
A primer design algorithm for Multiplex Polymerase Chain Reaction (PCR) and gene cloning. The algorithm works by first locating the first 30 base pairs downstream of the start codon for each gene in the reference genome (in this case Nematocida parisii). Any genes with such sequences below the optimal TM are filtered out. Otherwise, 1 bp is trimmed at a time from the end of the 30 bp sequence while the TM is above the optimal TM cutoff for PCR conditions. A similar process is done for the reverse strand. These are the initial DNA primer sequences. 

At this stage, a pool placement function sorts primers into subpools based on constraints including gene size, location in the genome, complementarity with primers currently in the pool, and the current size of the pool. If the primer is rejected, the next subpool is checked.

For all primers that pass these steps, a gateway, universal, and subpool-specific sequences are added on along with restriction sites to create the final oligonucleotides ready for PCR.

The final result is a pool of given number of subpools with equal size that contains oligos that are optimally designed to achieve high rates of success with PCR, and gene cloning.

The main module called PrimerDesign.py can be run directly to output an excel file of all successful candidate genes in the genome along with their oligo sequences. Pool size and number of subpools can be adjusted, as well as reference genome file, TM cutoff, and restriction sites.

## Making changes
