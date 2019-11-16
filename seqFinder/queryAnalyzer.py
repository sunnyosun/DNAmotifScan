import os
import csv
from collections import namedtuple

from .seqFromGenome import SeqFromGenomeArgs

class QueryAnalyzer:

    def __init__(self, query_file):
        self.query_file = query_file
        self.Query = namedtuple(
            "Query", ["index", "gene_id", "chr", "tss", "strand"])
    
    def get_queries(self):
        csv_f = open(self.query_file, "r")
        csv_f.readline()
        return map(self.Query._make, csv.reader(csv_f))


class QueryAnalyzerArgs(SeqFromGenomeArgs):

    def __init__(self):
        super().__init__()
        self.add_argument("-q", "--query-file", required=True,
                          help="Filename of the query file.")
