import os
import argparse
import pickle as pk

from Bio import SeqIO
from tqdm import tqdm

from .utils import get_motifs
from .findMotif import MotifFinder

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--motif-path",
                        help="Path to the directory storing motif files.")
    parser.add_argument("-f", "--fasta-path",
                        help="Path to the fasta file of sequences.")
    args = parser.parse_args()

    motifs = get_motifs(args.motif_path)
    finder = MotifFinder(motifs)
    sequences = SeqIO.parse(open(args.fasta_path, "r"), "fasta")

    print("Searching for matches...")
    matches = dict()
    counter = 0
    for seq in tqdm(sequences, total=83051):
        matches[seq.id] = finder.get_matching_motifs(seq)
        counter += 1
        if counter == 5:
            break
    
    key = list(matches.keys())[0]
    print(matches[key])
    os.makedirs("outputs", exist_ok=True)
    with open("outputs/matches", "wb") as f:
        pk.dump(matches, f)

    # df = optimize_matches(matches)
    # df.to_csv("outputs/matches.csv")
    