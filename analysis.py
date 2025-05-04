import argparse
import os
import numpy as np
from typing import List, Union
import symnmf as symnmf_py
from kmeans import kmeans
import sklearn


def parse() -> argparse.Namespace:
    """Parsing Function for user input

    Returns:
        argparse.Namespace: Parsing object to pass user input to main program
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('K', type=str)
    parser.add_argument('file_name', type=str)
    return parser.parse_args()


def read_file(filepath: str) -> List[List[float]]:
    """Function to read datapoints from file

    Args:
        filepath (str): filepath to .txt file containing valid datapoints

    Returns:
        List[List[float]: 2D Array, each element is a datapoint which is itself an array of floats
    """
    points: List[List[float]] = []
    with open(filepath, "r", encoding="utf-8") as file:
        for line in file:
            line = line.strip()
            point = line.split(",")
            try:
                point = list(map(float, point))
            except ValueError:
                print("An Error Has Occurred")
                exit(1)
            points.append(point)
    return points


def kmeans_silhouette_score(K: int, datapoints: List[List[float]], filepath: str) -> float:
    """Calculates KMeans (implementation from original HW1) silhouette score

    Args:
        K (int): Number of clusters
        datapoints (List[List[float]]): 2D array of datapoints 
        filepath (str): filepath to .txt file containing datapoints 

    Returns:
        float: KMeans silhouette score
    """
    kmeans_clusters = kmeans(K, 300, filepath)
    kmeans_pairwise_distances = sklearn.metrics.pairwise_distances(datapoints, kmeans_clusters)
    labels = np.argmin(kmeans_pairwise_distances, axis=1)
    silhouette_score = sklearn.metrics.silhouette_score(datapoints, labels)
    return silhouette_score
    

def symnmf_silhouette_score(K, datapoints):
    """Calculates SymNMF silhouette score

    Args:
        K (int): Number of clusters
        datapoints (List[List[float]]): 2D array of datapoints 

    Returns:
        float: SymNMF silhouette score
    """
    H = symnmf_py.nmf(K, datapoints)
    labels = np.argmax(H, axis=1)
    if len(np.unique(labels)) == 1: # Exception is raised if all datapoints assigned to same cluster
        print("An Error Has Occurred")
        exit()
    silhouette_score = sklearn.metrics.silhouette_score(datapoints, labels)
    return silhouette_score


def main():
    """Main function for analysis.py, shows output of silhouette comparison
    """
    args: argparse.Namespace = parse()
    K, file_name = args.K, args.file_name
    try:
        K = float(K)
        if K != int(K):
            raise ValueError()
        K = int(K)
    except ValueError:
        print("An Error Has Occurred")
        return
    filepath = os.path.join(os.path.join(os.getcwd()), file_name)

    points = read_file(filepath)
    if K <= 1 or K >= len(points):
        print("An Error Has Occurred")
        return
    
    kmeans_score = kmeans_silhouette_score(K, points, filepath)
    nmf_score = symnmf_silhouette_score(K, points)
    print(f"nmf: {nmf_score:.4f}")
    print(f"kmeans: {kmeans_score:.4f}")
    

if __name__ == "__main__":
    main()