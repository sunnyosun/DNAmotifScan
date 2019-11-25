import os
from collections import namedtuple
import statistics as st
import pickle as pk

import pandas as pd
import numpy as np
from tqdm import tqdm


class ResultDict(dict):

    def statistics(self):
        stat = dict()
        for seq_name, motifs in tqdm(self.items()):
            scores = list()
            for motif in motifs.values():
                scores += motif["scores"]

            stat[seq_name] = dict(
                score_mean=st.mean(scores),
                score_stdev=st.stdev(scores),
                n_match=len(scores)
            )
        return stat


class ResultParser:

    @classmethod
    def _combine(cls, motifs: list):
        buffer = dict()
        for motif in motifs:
            tf_name = "_".join(motif["motif"].split("_")[:-1])
            try:
                buffer[tf_name]["scores"] += motif["scores"]
                buffer[tf_name]["indices"] += motif["indices"]
                buffer[tf_name]["frequency"] += motif["frequency"]
            except KeyError:
                motif.pop("motif")
                buffer[tf_name] = motif
        return buffer

    @classmethod
    def _combine2(cls, motifs: dict):
        buffer = dict()
        for motif_name, data in motifs.items():
            tf_name = "_".join(motif_name.split("_")[:-1])
            try:
                buffer[tf_name]["scores"] += data["scores"]
                buffer[tf_name]["indices"] += data["indices"]
                buffer[tf_name]["frequency"] += data["frequency"]
            except KeyError:
                buffer[tf_name] = data
        return buffer

    @classmethod
    def _all_motifs(cls, loaded_dict):
        motifs = set()
        for data in loaded_dict.values():
            motifs = motifs.union(set(data.keys()))
        return list(motifs)

    @classmethod
    def _rearrange(cls, data):
        rearranged = ResultDict()
        for seq_name, motifs in data.items():
            motifs_dict = dict()
            for motif in motifs:
                motif_name = motif["motif"]
                motif.pop("motif")
                motifs_dict[motif_name] = motif
            rearranged[seq_name] = motifs_dict
        return rearranged

    @classmethod
    def load_results(cls, path, combine=True):
        results = list(os.scandir(path))
        raw_data = ResultDict()
        for result in tqdm(results):
            with open(result.path, "rb") as f:
                data = pk.load(f)
            if combine:
                raw_data.update(cls.combine_motifs(data))
            else:
                raw_data.update(cls._rearrange(data))
        return raw_data

    @classmethod
    def combine_motifs(cls, raw_data):
        combined_dict = ResultDict()
        for seq_name, motifs in raw_data.items():
            try:
                info = cls._combine(motifs)
            except Exception:
                info = cls._combine2(motifs)
            combined_dict[seq_name] = info
        return combined_dict

    @classmethod
    def frequency_matrix(cls,
                         loaded_dict=None,
                         combine_motifs=True,
                         columns=None):
        if combine_motifs:
            print("Combining motifs...")
            try:
                loaded_dict = cls.combine_motifs(loaded_dict)
            except KeyError:
                print("Motifs are combined already. Skipped combining.")
                pass
        if columns is None:
            columns = cls._all_motifs(loaded_dict)
        indices = list(loaded_dict.keys())

        row_dict = {idx:num for idx, num in zip(indices, range(len(indices)))}
        column_dict = {
            idx:num for idx, num in zip(columns, range(len(columns)))}
        contents = np.zeros((len(indices), len(columns))).tolist()
        print("Generating matrix...")
        for seq_name, motifs in tqdm(loaded_dict.items()):
            for motif, values in motifs.items():
                contents[row_dict[seq_name]][column_dict[motif]] = \
                    int(values["frequency"])
        frqc_mat = pd.DataFrame(
            contents, index=indices, columns=columns, dtype=int)
        return frqc_mat
