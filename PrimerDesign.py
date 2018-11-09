"""
primer script
"""
from Bio import SeqIO
from Bio.SeqUtils import Seq
import random
import time
from Primer import Primer


class Gene:
    """
    Structure for all the info for each gene.
    """

    def __init__(self, seq, start, end, length) -> None:
        """
        Initialize gene
        """
        self.start = int(start)
        self.end = int(end)
        self.length = int(length)
        self.seq = str(seq)

    def find_start(self) -> str:
        """
        Locates position of start codon, returns string of 30bp downstream from it.
        If no start codon, returns "no ATG"
        """
        try:
            start = self.seq.index("ATG")
            return self.seq[start:start+30]
        except ValueError:
            print("no ATG")

    def find_end(self) -> str:
        """
        last 30 bp of coding sequence in reverse complement
        """
        r = self.seq[::-1]
        stop_codons = ["GAT", "AAT", "AGT"]
        end = 1000000
        for codon in stop_codons:
            try:
                stop = r.index(codon)
                if stop < end:
                    end = stop
            except ValueError:
                end = -1
        if end >= 0:
            end += 3
        else:
            end = 0
        return str(Seq(r[end:end+30]).complement())


def get_distance(gene1: Gene, gene2: Gene, position: dict) -> int:
    """
    Returns distance between 2 genes in genome
    """
    return abs(position[gene1] - position[gene2])


def is_complementary(primer1, primer2) -> bool:
    """
    check for complementary of 10 bp or more
    """
    # some function to get reverse complement primer2

    primer2 = str(Seq(primer2).reverse_complement())
    start1, end1 = 0, 12
    while end1 <= len(primer1):
        start2, end2 = 0, 12
        while end2 <= len(primer2):
            if primer1[start1:end1] == primer2[start2:end2]:
                return True
            else:
                start2 += 1
                end2 += 1
        start1 += 1
        end1 += 1
    return False


def get_gene_length_diff(gene1: Gene, gene2: Gene) -> int:
    """
    Return diff in gene length
    """
    return gene1.length - gene2.length


def is_acceptable(query_gene, subpool, primer_set, gene_set, position: dict) -> bool:
    """
    Returns true if all the conditions of a primer against a pool are satisfied
    """

    for gene in subpool:
        if not(get_distance(query_gene,
                            gene, position) >= 10000 and get_gene_length_diff
               (gene_set[query_gene], gene_set[gene]) <= 500):
            return False
        if is_complementary(primer_set[query_gene][0], primer_set[gene][0]) \
                or is_complementary(primer_set[query_gene][0], primer_set[gene][1]) \
                or is_complementary(primer_set[query_gene][1], primer_set[gene][0]) \
                or is_complementary(primer_set[query_gene][1], primer_set[gene][1]):
            return False
    return True


def pool_placement(pool_size: int, tm_cutoff: int):
    """
    running algorithm to place primers in subpools
    """
    # gene_seqs_file = input("Enter gene sequences filename: ")
    # gene_table_file = input("Enter gene table file: ")
    # pool_size = int(input("Enter number of subpools desired: "))
    # Part 1 #
    # For each gene in the list we have a dictionary of key (gene ID) and tuple
    # containing gene start, end, length, first 30 bp, last 30 bp
    # Part 2 #

    # We check if any of the fwd or reverse primers from the gene dictionary
    # have 8 bp complementary. If any, we ensure they cannot be in same subpool.
    # Also we have to create the pool with 25 subs each containing 100 primers.
    fasta_sequences = \
        SeqIO.parse(open(file="GCF_000250985.1_Nema_parisii_ERTm1_V3_rna.fna"),
                    'fasta')
    geninfo, good_primers = {}, {}
    prot_file = open(file="ProteinTableTrim.txt")
    lines = prot_file.readlines()
    lines = lines[1:]
    gene_table = []
    for line in lines:
        gene_table.append(line.strip().split("\t"))
    i = 0
    for fasta in fasta_sequences:
        if i == len(gene_table):
            break
        # start, end, length are columns from the protein table
        geninfo[gene_table[i][4]] = Gene(str(fasta.seq), gene_table[i][1], gene_table[i][2], gene_table[i][5])
        good_primers[gene_table[i][4]] = []
        i += 1
    # Make a database to keep track of relative positions of all genes in the genome
    position = {}
    for i in range(len(gene_table) - 1):
        accession = gene_table[i][0]
        next_accession = gene_table[i + 1][0]
        curr_gene = gene_table[i][4]
        next_gene = gene_table[i + 1][4]
        if len(position) == 0:
            position[curr_gene] = int(gene_table[0][1])
        if accession != next_accession:
            position[next_gene] = position[curr_gene] + geninfo[curr_gene].end \
                                  - geninfo[curr_gene].start + \
                                  int(gene_table[i + 1][1])
        else:
            position[next_gene] = position[curr_gene] \
                                  + (int(int(gene_table[i + 1][1]) - int(gene_table[i][2]))) \
                                  + geninfo[curr_gene].end - geninfo[curr_gene].start

    # Algorithm for finding optimal fwd & rev primer based on TM#
    failed_genes = []
    for gene in geninfo:
        fwd = Primer(gene, "F", geninfo[gene].find_start())
        rev = Primer(gene, "R", geninfo[gene].find_end())
        if fwd.calc_tm() < tm_cutoff or rev.calc_tm() < tm_cutoff:
            good_primers.pop(gene)
            failed_genes.append(gene)
        else:
            while len(fwd.seq) > 18 and fwd.calc_tm() > tm_cutoff:
                fwd.seq = fwd.seq[:-1]
            good_primers[gene].append(fwd)
            while len(fwd.seq) > 18 and rev.calc_tm() > tm_cutoff:
                rev.seq = rev.seq[:-1]
            good_primers[gene].append(rev)

    # Brute force algorithm for placing primers into subpools #
    # 1. If a pool is empty, just add the primer
    # 2. If a pool already contains some primer(s):
    # a) check occupancy (<100?) of pool
    # b) check lengths (+/- 500?) against each primer in pool
    # c) check distance (>= 10kb?) against each primer in pool
    # d) check complementarity (<= 8bp) against each primer in pool
    # 4. If ALL conditions satisfied, add primer to the pool, else skip pool
    # 5. If this is the last pool and fails, add to failed genes list
    print(len(failed_genes))
    keys_list = list(good_primers.keys()).copy()
    keys_list = keys_list[50:150]
    num_genes = len(keys_list)
    # Make big pool #
    num_pools = num_genes // pool_size + 1
    super_pool = [{} for _ in range(num_pools)]
    while len(keys_list) != 0:
        if len(keys_list) == 1:
            print("lol")
        curr_gene = random.choice(keys_list)
        fwd = good_primers[curr_gene][0]
        rev = good_primers[curr_gene][1]
        i = 0
        while i < num_pools:
            if is_complementary(fwd, rev):
                failed_genes.append(curr_gene)
                keys_list.remove(curr_gene)
                break
            curr_pool = super_pool[i]
            if len(curr_pool) == 0 or \
                    (is_acceptable(curr_gene, curr_pool, good_primers, geninfo, position) and
                     len(curr_pool) < pool_size):
                curr_pool[curr_gene] = (fwd, rev)
                keys_list.remove(curr_gene)
                break
            else:
                i += 1
                if i == num_pools:
                    failed_genes.append(curr_gene)
                    keys_list.remove(curr_gene)
    print(len(failed_genes))
    # insert (in front) gateway f/r, universal f, subpool f.
    # append Sap/Bsa, subpool r.
    return super_pool, good_primers


def add_subpoolf(oligo, super_pool, kmers):
    """
    # do we still include the universal subpool f for size 25 or 50
    # do we want to have the same subpool
    """
    num = 1
    i = 0
    while i < num:
        kmer = random.choice(kmers)
        for gene in super_pool:
            if is_complementary(super_pool[gene][0], kmer) or \
               is_complementary(super_pool[gene][1], kmer):
                kmer = random.choice(kmers)
        oligo[0].insert(0, kmer)
        oligo[1].insert(0, kmer)
        kmers.remove(kmer)
        i += 1
    return oligo


def add_universalf(oligo, super_pool, kmers):
    """
    stuff
    """
    # randomly select a 20-mer
    kmer = random.choice(kmers)
    for gene in super_pool:
        if is_complementary(super_pool[gene][0], kmer) or \
           is_complementary(super_pool[gene][1], kmer):
            kmer = random.choice(kmers)
    oligo[0].insert(0, kmer)
    oligo[1].insert(0, kmer)
    kmers.remove(kmer)
    return oligo


def add_gateway(oligo):
    """
    adds gateway seq
    """
    gatewayf = "ACAAGTTTGTACAAAAAAGCAGGCTCA"
    gatewayr = "ACCACTTTGTACAAGAAAGCTGGGTT"

    # check complementary gateway f and r with primers

    oligo[0].insert(0, gatewayf)
    oligo[1].insert(0, gatewayr)
    return oligo


def add_restriction_sites(oligo):
    """
    screen & append restriction site after the specific gene sequence
    """
    sap_seq = "GGTCGAAGAGC"
    bsa_seq = "ACTGAGAGACC"
    sap1 = ("GCTCTTC", "GAAGAGC")
    bsa1 = ("GGTCTC", "GAGACC")
    if sap1[0] in oligo[0] or sap1[1] in oligo[0]:
        if bsa1[0] in oligo[0] or bsa1[1] in oligo[0]:
            oligo[0].append("no RS")
        else:
            oligo[0].append(bsa_seq)
    else:
        oligo[0].append(sap_seq)

    if sap1[0] in oligo[1] or sap1[1] in oligo[1]:
        if bsa1[0] in oligo[1] or bsa1[1] in oligo[1]:
            oligo[1].append("no RS")
        else:
            oligo[1].append(bsa_seq)
    else:
        oligo[1].append(sap_seq)
    return oligo


def add_subpoolr(oligo, super_pool, kmers):
    """
    screen & append subpool specific reverse from list of 20-mers to end
    """
    kmer = ""
    for gene in super_pool:
        kmer = random.choice(kmers[(3*len(kmers)//4):])
        if is_complementary(super_pool[gene][0], kmer) or \
                is_complementary(super_pool[gene][1], kmer):
            kmer = random.choice(kmers)
    oligo[0].append(kmer)
    oligo[1].append(kmer)
    kmers.remove(kmer)
    return oligo


def create_oligos(super_pool: list, rand_dna: list):
    """
    Create good_oligos by adding the subpool, universal etc.
    """

    good_oligos = {}
    for subpool in super_pool:
        for gene in subpool:
            good_oligos[gene] = [[subpool[gene][0]], [subpool[gene][1]]]

    for subpool in super_pool:
        for gene in subpool:
            good_oligos[gene] = add_gateway(good_oligos[gene])
            good_oligos[gene] = add_universalf(good_oligos[gene], subpool, rand_dna.copy())
            good_oligos[gene] = add_subpoolf(good_oligos[gene], subpool, rand_dna.copy())
            good_oligos[gene] = add_restriction_sites(good_oligos[gene])
            good_oligos[gene] = add_subpoolr(good_oligos[gene], subpool, rand_dna.copy())
    print("good_oligos done")
    return good_oligos


if __name__ == '__main__':
    t0 = time.time()

    kmer_file = "DNA_20mers.txt"
    f = open(kmer_file)
    flines, kmer_list = f.readlines(), []
    for x in range(1, len(flines), 2):
        kmer_list.append(flines[x].strip().split("\t")[0])

    pool, primers = pool_placement(100, 50)

    oligos = create_oligos(pool, kmer_list)

    for oligonucleotide in oligos:
        print(oligonucleotide)
    print(time.time() - t0)


# outer distances (first part of 1st gene ending of 2nd gene)
# 8 base window check instead of blast (since blast includes gaps)
# see how primer programs have solved this problem
# ungapped_alignment = yes option to blast without gaps
# http://biopython.org/DIST/docs/api/Bio.pairwise2-module.html customize gapping


# 21% of the primers ruled out with TM cutoff 55
# an additional 6% have at least 8bp complementary w/something
# program run time: 3:30 - 4:00 minutes
# do all the pool sizes, and change the # of bp compl. and the TM
# check each constraint independently


# screen sap1 sites
# 5' G C T C T T C 3' design universal sequence (TM >= 55, GC content, 200 pool size)

# % that don't map to the genome blast them against some db (uniprot)
# self map all THE NEW GENOMES
#


# check the contamination  of C elegans and E coli for the new speices
# Epiphaga, Minor, not grown in C Elegans
#
