import os
from collections import namedtuple
import statistics as st

Info = namedtuple(
        "Info",
        "motifs best_match n_motifs score_avg score_stdev"
    )


class ResultParser:

    @classmethod
    def _unwrap(cls, buffer):
        motifs = list()
        scores = list()
        best_score = 0
        best_match = list()
        for motif_name, data in buffer.items():
            motifs.append((motif_name, data["score"], data["index"]))
            scores.append(data["score"])
            if data["score"] > best_score:
                best_score = data["score"]
                best_match = list()
                best_match.append((motif_name, data["score"], data["index"]))
            elif data["score"] == best_score:
                best_match.append((motif_name, data["score"], data["index"]))
        n_motifs = len(motifs)
        score_avg = st.mean(scores)
        score_stdev = st.stdev(scores)
        return Info(motifs, best_match, n_motifs, score_avg, score_stdev)

    @classmethod
    def _refine(cls, motifs: list):
        buffer = dict()
        for motif in motifs:
            motif_name = "_".join(motif["motif"].split("_")[:-1])
            if motif_name not in buffer:
                buffer[motif_name] = motif
            else:
                if buffer[motif_name]["score"] < motif["score"]:
                    buffer[motif_name] = motif
        return cls._unwrap(buffer)

    @classmethod
    def parse(cls, results):
        parsed_dict = dict()
        for seq_name, motifs in results.items():
            info = cls._refine(motifs)
            parsed_dict[seq_name] = info
        return parsed_dict
