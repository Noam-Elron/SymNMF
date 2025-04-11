import argparse
import os
import numpy as np
import math
from typing import List, Union
import symnmf as symnmf_py
from kmeans_hw1 import kmeans
import sklearn


def parse() -> argparse.Namespace:
    """

    Parameters
    ----------
    distances_list : _type_
        _description_
    total_distances : _type_
        _description_
        
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('K', type=str)
    parser.add_argument('file_name', type=str)
    return parser.parse_args()


def read_file(filepath):
    points = []
    with open(filepath, "r", encoding="utf-8") as file:
        for line in file:
            line = line.strip()
            point = line.split(",")
            point = list(map(float, point))
            points.append(point)
    return points


def kmeans_silhouette_score(K, datapoints, filepath):
    kmeans_clusters = kmeans(K, 300, filepath)
    kmeans_pairwise_distances = sklearn.metrics.pairwise_distances(datapoints, kmeans_clusters)
    labels = np.argmin(kmeans_pairwise_distances, axis=1)
    silhouette_score = sklearn.metrics.silhouette_score(datapoints, labels)
    return silhouette_score
    

def symnmf_silhouette_score(K, datapoints):
    H = symnmf_py.nmf(K, datapoints)
    labels = np.argmax(H, axis=1)
    silhouette_score = sklearn.metrics.silhouette_score(datapoints, labels)
    return silhouette_score


def main():
    args: argparse.Namespace = parse()
    K, file_name = args.K, args.file_name
    K = int(K)
    filepath = os.path.join(os.path.join(os.getcwd()), file_name)

    points = read_file(filepath)
    kmeans_score = kmeans_silhouette_score(K, points, filepath)
    nmf_score = symnmf_silhouette_score(K, points)
    print(f"nmf: {nmf_score:.4f}")
    print(f"kmeans: {kmeans_score:.4f}")
    

if __name__ == "__main__":
    main()