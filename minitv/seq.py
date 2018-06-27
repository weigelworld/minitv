"""
Handles sequence file processing
"""

from Bio import SeqIO


class FASTAFile:
    """
    Represents a FASTA file of one or more sequences.
    """

    def __init__(self, path):
        """

        :param path: path to the FASTA file
        :type path: pathlib.Path
        """

        self.path = path
        self.name = path.stem

        self.records = list(SeqIO.parse(str(path), 'fasta'))

    def __len__(self):
        """

        :return: number of records in the FASTA file
        """

        return len(self.records)
