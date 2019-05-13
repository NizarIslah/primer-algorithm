"""
Representation of a gene.
"""
from Bio.Seq import Seq


class Gene:
    """
    Structure for all the info for each gene.
    """

    def __init__(self, seq, start, end, length, name) -> None:
        """
        Initialize gene
        """
        self.start = int(start)
        self.end = int(end)
        self.length = int(length)
        self.seq = str(seq)
        self.name = name

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
        TODO: This should instead search from the first ATG until a stop codon is reached.
        """
        rev = self.seq[::-1]
        stop_codons = ["GAT", "AAT", "AGT"]
        end = 1000000
        for codon in stop_codons:
            try:
                stop = rev.index(codon)
                if stop < end:
                    end = stop
            except ValueError:
                end = -1
        if end >= 0:
            end += 3
        else:
            end = 0
        return str(Seq(rev[end:end+30]).complement())

    def get_length_diff(self, gene) -> int:
        """
        Return diff in gene length
        """
        if isinstance(gene, Gene):
            return self.length - gene.length

    def get_distance(self, gene, position: dict) -> int:
        """
        Returns distance between 2 genes in genome.
        """
        return abs(position[self.name] - position[gene.name])

    def __str__(self) -> str:
        """
        String representation of a gene
        """
        return "Name: {}, Length: {}".format(self.name, self.length)
