# primer-algorithm

## Overview
A primer design algorithm for Multiplex Polymerase Chain Reaction (PCR) and gene cloning. The algorithm works by first locating the first 30 base pairs downstream of the start codon for each gene in the reference genome (in this case Nematocida parisii). Any genes with such sequences below the optimal TM are filtered out. Otherwise, 1 bp is trimmed at a time from the end of the 30 bp sequence while the TM is above the optimal TM cutoff for PCR conditions. A similar process is done for the reverse strand.

At this stage, a pool placement function sorts primers into subpools based on constraints including gene size, location in the genome, complementarity with primers currently in the pool, and the current size of the pool. If the primer is rejected, the next subpool is checked.

For all primers that pass these steps, a gateway, universal, and subpool-specific sequences are added on along with restriction sites to create the final oligonucleotides ready for initial round of PCR.

The final result is a batch of equally sized subpools containing oligos that are optimally designed to achieve high rates of success with PCR, and gene cloning.

## Running the program

The main module called Main.py can be run directly to output an excel file of all potential candidate genes in the genome along with their designed oligo sequences. Pool size and number of subpools can be adjusted, as well as reference genome file, TM cutoff, and restriction sites. These changes can all be made in the Main.py file.
