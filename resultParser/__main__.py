import os
import pickle as pk
from .resultParser import ResultParser

from tqdm import tqdm


if __name__ == "__main__":
    parsed = dict()
    outputs = list(os.scandir("temp/outputs"))
    for output in tqdm(outputs):
        with open(output.path, "rb") as f:
            results = pk.load(f)
        parsed.update(ResultParser.parse(results))
    
    with open("data/parsed_results.pk", "wb") as f:
        pk.dump(parsed, f)


    # with open("temp/outputs/matches_splitted_0", "rb") as f:
    #     results = pk.load(f)
    # # rp = ResultParser()
    # parsed = ResultParser.parse(results)
    # key1 = list(parsed.keys())[0]
    # print("sequence name: {}".format(key1))
    # info = parsed[key1]
    # print("motifs: {}".format(info.motifs))
    # print("best match: {}".format(info.best_match))
    # print("n_motifs: {}".format(info.n_motifs))
    # print("score_avg: {}".format(info.score_avg))
    # print("score_stdev: {}".format(info.score_stdev))
