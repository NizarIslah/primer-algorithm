"""
primer class
"""
from Bio.SeqUtils import MeltingTemp as Mt, Seq


class Primer:
    """
    class for primer
    """
    def __init__(self, gene_id: str, direction: str, seq: str, tm=0) -> None:
        """
        initialize a primer
        """
        self.gene_id = gene_id
        self.direction = direction
        self.seq = seq
        self.tm = tm

    def calc_tm(self) -> float:
        """
        Calculates the TM of the primer sequence.
        """
        self.tm = float('%0.2f' % Mt.Tm_NN(Seq(self.seq), nn_table=Mt.DNA_NN2))
        return self.tm

    def __str__(self) -> str:
        """
        str rep
        """
        return "ID: {} \n Direction: {} \n Sequence: {}".format(self.gene_id, self.direction, self.seq)
