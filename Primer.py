"""
primer class
"""


class Primer:
    """
    class for primer
    """
    def __init__(self, direction: str, template: str, seq: str) -> None:
        """
        initialize a primer
        """
        self.direction = direction
        self.template = template
        self.seq = seq

    def __str__(self) -> str:
        """
        str rep
        """
        return "direction: {}, seq: {}".format(self.direction, self.seq)
