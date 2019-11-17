import os

import numpy as np

def get_motifs(path):
    files = os.scandir(path)
    motifs = dict()
    for f in files:
        if not f.is_file():
            continue
        motifs[f.name] = np.fromfile(f.path, sep="\t").reshape(-1, 4)
    return motifs
