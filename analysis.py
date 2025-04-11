import argparse
import os
import numpy as np
import math
from typing import List, Union
import symnmf
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

def main():
    args: argparse.Namespace = parse()
    K, file_name = args.K, args.file_name
    
    try:
        K = int(K)
        filepath = os.path.join(os.path.join(os.getcwd()), file_name)

    except ValueError:
        print("An Error Has Occurred py1")
    

    points = read_file(filepath)

    if K <= 1 or K >= len(points):
        print("An Error Has Occurred py2")
        return
    
    kmeans_clusters = kmeans(K, 600, filepath)
    kmeans_pairwise_distances = sklearn.metrics.pairwise_distances(points, kmeans_clusters)
    
    kmeans_silloute_score = sklearn.metrics.silhouette_score()

