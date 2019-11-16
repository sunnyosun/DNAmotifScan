""" Find sequence from genome fasta file by position
"""

import os
import gzip
from functools import partial
import logging
from argparse import ArgumentParser

from Bio import SeqIO


class FileTypeError(Exception):
    pass

class SeqFromGenome:

    def __init__(self, genome_file="hg38.fa.gz"):
        self.genome_f = os.path.join("data", genome_file)
        self.gz = self._check_file_type(genome_file)
        self.genome = self._read_genome()
        self.logger = self._getLogger()

    def _check_file_type(self, file_name):
        if file_name.endswith(".gz"):
            return True
        elif file_name.endswith(".fa") or file_name.endswith(".fasta"):
            return False
        else:
            raise FileTypeError("Can not infer the filetype automatically. "
                "Please use .fa/.fasta/.gz for the file extension.")

    def _read_genome(self):
        if self.gz:
            _open = partial(gzip.open, mode="rt")
        else:
            _open = open
        with _open(self.genome_f) as f:
            contents = SeqIO.to_dict(SeqIO.parse(f, "fasta"))
        return contents

    def find_sequence(self,
                      chromosome,
                      position,
                      strand,
                      upperstream=1000,
                      downstream=200):

        if strand == "-":
            upperstream, downstream = downstream, upperstream
        seq = self.genome[chromosome][
            position-upperstream:position+downstream+1].upper()
        N_occu = seq.seq.count("N") / float(len(seq))
        if N_occu > 0.5:
            self.logger.warning("Chromosome: {}\nPosition: {}\nStrand: "
                "{}\nUpperstream: {}\nDownstream: {}\nThe seqence has more "
                "than 50% Ns.".format(
                        chromosome, position, strand, upperstream, downstream))
        return seq

    def _getLogger(self):
        logger = logging.getLogger(self.__class__.__name__)
        hl = logging.StreamHandler()
        formmatter = logging.Formatter("%(asctime)-15s %(message)s")
        hl.setFormatter(formmatter)
        logger.addHandler(hl)
        return logger


class SeqFromGenomeArgs(ArgumentParser):

    def __init__(self):
        super().__init__()
        self.add_argument("-f", "--fasta-filename", default="hg38.fa.gz",
                          help="Filename of the fasta file.")
        self.add_argument("-u", "--upperstream", type=int, default=1000,
                          help="Upperstream basepairs.")
        self.add_argument("-d", "--downstream", type=int, default=200,
                          help="Downstream basepais.")
