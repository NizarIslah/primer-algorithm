"""
Main primer design module.
"""
import random
from typing import Union
from Bio import SeqIO
from Bio.SeqUtils import Seq
from Primer import Primer
from Gene import Gene
import matplotlib.pyplot as plt
import csv


def is_complementary(primer1: Union[Primer, str], primer2: Union[Primer, str], bp=12) -> bool:
    """
    Check for complementarity of bp diff or more between 2 primers or sequences.
    """
    p1seq = str(Seq(primer1.seq)) if isinstance(primer1, Primer) else primer1
    p2seq = str(Seq(primer2.seq).reverse_complement()) if isinstance(primer2, Primer) else primer2
    start1, end1 = 0, bp
    while end1 <= len(p1seq):
        start2, end2 = 0, bp
        while end2 <= len(p2seq):
            if p1seq[start1:end1] == p2seq[start2:end2]:
                return True
            else:
                start2 += 1
                end2 += 1
        start1 += 1
        end1 += 1
    return False


def is_acceptable(query_gene, subpool, primer_set, gene_set, position: dict) -> bool:
    """
    Returns true if all the conditions of a primer against a pool are satisfied
    """

    for gene in subpool:
        if not(gene_set[query_gene].get_distance(gene_set[gene], position) >= 10000 and
               (gene_set[query_gene].get_length_diff(gene_set[gene])) <= 500):
            return False
        if is_complementary(primer_set[query_gene][0], primer_set[gene][0]) \
                or is_complementary(primer_set[query_gene][0], primer_set[gene][1]) \
                or is_complementary(primer_set[query_gene][1], primer_set[gene][0]) \
                or is_complementary(primer_set[query_gene][1], primer_set[gene][1]):
            return False
    return True


def filter_tm(good_primers: dict, geninfo: dict, tm_cutoff: int) -> (dict, list):
    """
    Filters out primers below the tm cutoff
    """
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
    return good_primers, failed_genes


def load_genes(num=2660, seqfile="GCF_000250985.1_Nema_parisii_ERTm1_V3_rna.fna", protfile="ProteinTableTrim.txt"):
    """
    Create gene dictionary
    """
    fasta_sequences = SeqIO.parse(open(file=seqfile), 'fasta')
    geninfo, init_primers = {}, {}
    prot_file = open(file=protfile)
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
        geninfo[gene_table[i][4]] = Gene(str(fasta.seq), gene_table[i][1], gene_table[i][2],
                                         gene_table[i][5], gene_table[i][4])
        init_primers[gene_table[i][4]] = []
        i += 1
    # Make a dictionary to keep track of relative positions of all genes in the genome
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
    subset_primers = {g: init_primers[g] for g in list(init_primers.keys())[:num]}
    subset_geninfo = {g: geninfo[g] for g in list(geninfo.keys())[:num]}
    return subset_primers, subset_geninfo, position


def pool_placement(pool_size: int, tm_cutoff: int, init_primers, geninfo, position):
    """ Place into pools. """
    # Algorithm for finding optimal fwd & rev primer
    # 1. If a pool is empty, just add the primer
    # 2. If a pool already contains some primer(s):
    # a) If ALL conditions satisfied, add primer to the pool, else skip pool
    # b) If fails to meet conditions for all subpools, add to failed genes list
    good_primers, failed_genes = filter_tm(init_primers, geninfo, tm_cutoff)
    keys_list = list(good_primers.keys()).copy()
    num_genes = len(keys_list)
    # success = 0
    num_pools = num_genes // pool_size + 1
    super_pool = [{} for _ in range(num_pools)]
    while len(keys_list) != 0:
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
                # success += 1
                keys_list.remove(curr_gene)
                break
            else:
                i += 1
                if i == num_pools:
                    failed_genes.append(curr_gene)
                    keys_list.remove(curr_gene)
    # success = success/num_genes*100
    return super_pool


def add_subpoolf(oligo, super_pool, kmers):
    """
    Insert the forward subpool specific sequence.
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
    Insert the universal forward sequence.
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
    Insert the gateway sequence.
    """
    gatewayf = "ACAAGTTTGTACAAAAAAGCAGGCTCA"
    gatewayr = "ACCACTTTGTACAAGAAAGCTGGGTT"

    # check complementary gateway f and r with primers

    oligo[0].insert(0, gatewayf)
    oligo[1].insert(0, gatewayr)
    return oligo


def add_restriction_sites(oligo):
    """
    Screen & append restriction site after the specific gene sequence.
    """
    sap_seq = "GGTCGAAGAGC"
    bsa_seq = "ACTGAGAGACC"
    sap1 = ("GCTCTTC", "GAAGAGC")
    bsa1 = ("GGTCTC", "GAGACC")
    if sap1[0] in oligo[0] or sap1[1] in oligo[0]:
        if bsa1[0] in oligo[0] or bsa1[1] in oligo[0]:
            oligo = "FAILED"
        else:
            oligo[0].append(bsa_seq)
    else:
        oligo[0].append(sap_seq)

    if sap1[0] in oligo[1] or sap1[1] in oligo[1]:
        if bsa1[0] in oligo[1] or bsa1[1] in oligo[1]:
            oligo = "FAILED"
        else:
            oligo[1].append(bsa_seq)
    else:
        oligo[1].append(sap_seq)
    return oligo


def add_subpoolr(oligo, super_pool, kmers):
    """
    Screen & append subpool specific reverse from list of 20-mers to end.
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


def list_from_file(filename):
    """ Make list of sequences from the file"""
    file = open(filename)
    flines, kmer_list = file.readlines(), []
    for x in range(1, len(flines), 2):
        kmer_list.append(flines[x].strip().split("\t")[0])
    return kmer_list


def create_oligos(super_pool: list, rand_dna_file):
    """
    Create good_oligos by adding the subpool, universal etc.
    """
    rand_dna = list_from_file(rand_dna_file)
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
    return good_oligos


def save_output(path, fmt, size, oligos):
    """ Save results to output file """
    out = path + str(size) + fmt
    with open(out, "w") as f:
        writer = csv.writer(f)
        for gene_id in oligos:
            result = get_oligo_seq(oligos, gene_id)
            writer.writerow([gene_id, result])
        f.close()


def get_oligo_seq(oligonucs: dict, geneid: str) -> str:
    """
    Prints the gene id and string representing the entire oligo sequence.
    """
    oligo_set = oligonucs[geneid]
    s = ''
    for oligo in oligo_set:
        for o in oligo:
            if isinstance(o, str):
                s += o
            else:
                s += o.seq
        s += "\t"
    return s


def make_line_plot(data, title, x_axis, y_axis, path):
    """ Makes a plot with the data."""
    fig = plt.figure()
    x, y = zip(*sorted(data.items(), key=lambda kv: kv[0]))
    plt.plot(x, y)
    fig.suptitle(title)
    plt.xlabel(x_axis)
    plt.ylabel(y_axis)
    fig.savefig(path)
    plt.show()
