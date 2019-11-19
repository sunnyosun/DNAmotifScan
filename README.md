# DNAmotifScan

This script scans through an input sequences and looks for possible TF binding sites based on motifs.
The motifs are provided as Position Weighted Matrix (pwm).
The sequece is provided in fasta format.

Examples are included.

### Create sequences file
`python -m seqFinder -q <path_to_csv_file> [-f path_to_fasta_file]`

### Find motifs
`python -m motifFinder -p path_to_motif_files_directory -f path_to_fastafile`

### Parse and save the found motifs
`python -m resultParser`

The save data is a dictionary of namedtuples. Keys are the names of the sequences. Namedtuple contains the parsed data, which has name fields: `motifs`, `best_match`, `n_motifs`, `score_avg`, `score_stdev`. `motifs` and `best_match` are list of tuples. Data in the tuples are motif name, score, index (position in the sequence)

The parsed data were saved with pickle. To correctly unpickle the data, the unpickle process must be done inside this package.
