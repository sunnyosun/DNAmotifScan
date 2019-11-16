from argparse import ArgumentParser

from .seqFromGenome import SeqFromGenome
from .queryAnalyzer import QueryAnalyzer, QueryAnalyzerArgs

if __name__ == "__main__":
    arg_parser = QueryAnalyzerArgs()
    arg_parser.add_argument("-o", "--output", default="sequences.fa",
                            help="Output file name. [default: sequences.fa]")
    args = arg_parser.parse_args()

    print("Reading fasta file...")
    seqFinder = SeqFromGenome(args.fasta_filename)
    qanalyzer = QueryAnalyzer(args.query_file)

    queries = qanalyzer.get_queries()
    output_f = open(args.output, "w")
    counter = 0
    for q in queries:
        if not q.chr.startswith("chr"):
            continue
        seq = seqFinder.find_sequence(
            q.chr, int(q.tss), q.strand,
            args.upperstream, args.downstream)
        output_f.write(
            ">"+\
            "_".join([q.chr, q.tss, q.strand, q.gene_id])+\
            "\n")
        output_f.write(str(seq.seq)+"\n")
        counter += 1
        if counter % 10 == 0:
            print(
                "{} sequences written into {}...".format(counter, args.output),
                end="\r")
    output_f.close()
