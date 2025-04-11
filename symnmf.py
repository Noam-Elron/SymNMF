import argparse
import os
import numpy as np
import math
from typing import List, Union
import symnmf


def pretty_print(matrix):
    for i in range(len(matrix)):
        print(",".join(map(lambda x:  f'{x:.4f}', matrix[i])))



def euclidean_distance(point, other) -> float:
    """

    Parameters
    ----------
    distances_list : _type_
        _description_
    total_distances : _type_
        _description_
        
    """
    total = 0
    for i in range(len(point)):
        total += pow((point[i] - other[i]), 2)
    return math.sqrt(total)

def initialize_H(norm_matrix, K):
    np.random.seed(1234)
    m = np.mean(norm_matrix)
    dimension = len(norm_matrix)
    H = []

    for i in range(dimension*K):
        val = np.random.uniform(0, 2 * math.sqrt(m / K))
        H.append(val)
    H = np.reshape(H, (dimension, K))
    return H.tolist()
    


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
    parser.add_argument('goal', type=str)
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

def sym(points):
    return symnmf.sym(points)

def diag(points):
    return symnmf.diag(points)

def norm(points):
    return symnmf.norm(points)

def nmf(K, points):
    norm_matrix = norm(points)
    H = initialize_H(norm_matrix, K)
    return symnmf.symnmf(H, norm_matrix)


def main():
    args: argparse.Namespace = parse()
    K, goal, file_name = args.K, args.goal, args.file_name
    
    try:
        K = int(K)
        filepath = os.path.join(os.path.join(os.getcwd()), file_name)

    except ValueError:
        print("An Error Has Occurred")
    

    points = read_file(filepath)

    if K <= 1 or K >= len(points):
        print("An Error Has Occurred")
        return

    goals_mapping = {"symnmf": 0, "sym": 1, "ddg": 2, "norm": 3}

    if goal not in goals_mapping:
        print("An Error Has Occurred")
        return

    if goals_mapping[goal] == 0:
        pretty_print(nmf(K, points))
    elif goals_mapping[goal] == 1:
        pretty_print(sym(points))
    elif goals_mapping[goal] == 2:
        pretty_print(diag(points))
    elif goals_mapping[goal] == 3:
        pretty_print(norm(points))



if __name__ == "__main__":
    main()

