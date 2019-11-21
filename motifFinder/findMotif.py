import numpy as np
import pandas as pd

class MotifFinder:

    def __init__(self, motifs: dict):
        self.pwms = {key:self._get_pwm(motif) for key, motif in motifs.items()}
        self.pwm_order = {'A':0,'C':1,'G':2,'T':3,'N':4}
        
    def _get_pwm(self, motif):
        # insert the 5th column for N
        return np.concatenate((motif, np.zeros((motif.shape[0], 1))), axis=1)

    def _seq2mat(self, sequence):
        mat = np.zeros((len(self.pwm_order), len(sequence)))
        for i, s in enumerate(sequence):
            mat[self.pwm_order[s], i] = 1
        return mat

    def _get_max_score(self, pwm):
        """ Get the score of the best matching sequence
        """
        return np.sum(pwm.max(axis=1))

    def _sum_diagnols(self, mat):
        sums = list()
        for i in range(mat.shape[1]-mat.shape[0]+1):
            sums.append(np.trace(mat[:, i:]))
        return np.array(sums).reshape(-1, 1)

    def _similarity_score(self, pwm, forward_mat, reverse_mat, threshold=0.9):
        """
        Calcualte the similarity score for a given sequence in each bp
        width of the sliding window == length of motif
        """
        best_match_score = self._get_max_score(pwm)

        forward_scores = self._sum_diagnols(np.dot(pwm, forward_mat))
        forward_scores = forward_scores / best_match_score
        reverse_scores = self._sum_diagnols(np.dot(pwm, reverse_mat))
        reverse_scores = np.flip(reverse_scores)
        reverse_scores = reverse_scores / best_match_score
        scores = np.concatenate((forward_scores, reverse_scores), axis=1)
        strands = (scores[:, 0] > scores[:, 1]).astype(int)
        best_scores = np.max(scores, axis=1)
        matched_indices = np.where(best_scores>threshold)
        matched_scores = best_scores[matched_indices].tolist()
        matched_indices = matched_indices[0].tolist()
        # strand = strands.squeeze()[best_index]
        frequency = np.sum(best_scores > threshold)
        return matched_scores, matched_indices, frequency

    def get_matching_motifs(self, sequence, threshold=0.9):
        matched_motifs = list()
        sequence = sequence.upper()
        forward_mat = self._seq2mat(sequence)
        reverse_mat = self._seq2mat(sequence.reverse_complement())
        for motif, pwm in self.pwms.items():
            scores, indices, frequency = self._similarity_score(
                                        pwm, forward_mat, reverse_mat)
            if len(scores) > 0:
                matched = dict()
                matched["motif"] = motif
                matched["scores"] = scores
                matched["indices"] = indices
                matched["frequency"] = frequency
                matched_motifs.append(matched)
        return matched_motifs


if __name__ == "__main__":
    from Bio import SeqIO
    from Bio.Alphabet import DNAAlphabet

    mf = MotifFinder("../data/pwm/AFP_1", "../sequences.fa")
    seq = next(mf.sequences).seq
    print("Length of sequence: {}".format(len(seq)))
    scores = mf.similarity_score(next(mf.sequences).seq)
    print("Scores: {}".format(scores))
    print("Shape of score array: {}".format(scores.shape))
    print("Shape of pwm matrix: {}".format(mf.pwm.shape))
