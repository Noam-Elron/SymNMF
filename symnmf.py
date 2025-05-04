import argparse
import os
import numpy as np
import math
from typing import List, Union
import symnmf_c
import sys


def pretty_print(matrix: List[List[float]]):
    """Function to print matrices as per project specifications 

    Args:
        matrix (List[List[float]]): Matrix we want to print.
    """
    for i in range(len(matrix)):
        print(",".join(map(lambda x:  f'{x:.4f}', matrix[i])))

def euclidean_distance(point: List[float], other: List[float]) -> float:
    """Calculates euclidean distance of two points

    Args:
        point (List[List[float]]): Point 1 we are trying to calculate euclidean distance to other point
        other (List[List[float]]): Point 2 we are trying to calculate euclidean distance to point

    Returns:
        float: euclidean distance of the two points
    """
    total = 0
    for i in range(len(point)):
        total += (point[i] - other[i]) * (point[i] - other[i])
    return math.sqrt(total)

def initialize_H(norm_matrix: List[List[float]], K: int) -> List[List[float]]:
    """Creates initial H matrix as per project specifications, uses np random seed 1234.

    Args:
        norm_matrix (List[List[float]]): Previously calculated norm matrix
        K (int): Number of clusters.

    Returns:
        List[List[float]]: Resultant initialized H matrix
    """
    np.random.seed(1234)
    m = np.mean(norm_matrix)
    dimension = len(norm_matrix)
    H = []
    upper_bound = 2 * np.sqrt(m / K)
    H = np.random.uniform(0, upper_bound, size=(dimension, K))
    return H.tolist()
    


def parse() -> argparse.Namespace:
    """Parsing Function for user input
    Returns:
        argparse.Namespace: Parsing object to pass user input to main program
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('K', type=str)
    parser.add_argument('goal', type=str)
    parser.add_argument('file_name', type=str)

    if len(sys.argv) != 4:
        print("An Error Has Occurred")
        exit(1)

    return parser.parse_args()


def read_file(filepath: str) -> List[List[float]]:
    """Function to read datapoints from file
    Args:
        filepath (str): filepath to .txt file containing valid datapoints
    Returns:
        List[List[float]]: 2D Array, each element is a datapoint which is itself an array of floats.
    """
    points = []
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

def sym(points: List[List[float]]) -> List[List[float]]:
    """Python wrapper function to calculate similarity matrix by calling appropriate C module function.
    Args:
        points (List[List[float]]): Datapoints used to calculate similarity matrix
    Returns:
        List[List[float]]: Resultant similarity matrix
    """
    return symnmf_c.sym(points)

def diag(points: List[List[float]]) -> List[List[float]]:
    """Python wrapper function to calculate diagonal matrix by calling appropriate C module function.
    Args:
        points (List[List[float]]): Datapoints used to calculate diagonal matrix
    Returns:
        List[List[float]]: Resultant diagonal matrix
    """
    return symnmf_c.diag(points)

def norm(points: List[List[float]]) -> List[List[float]]:
    """Python wrapper function to calculate norm matrix by calling appropriate C module function.
    Args:
        points (List[List[float]]): Datapoints used to calculate norm matrix
    Returns:
        List[List[float]]: Resultant norm matrix
    """
    return symnmf_c.norm(points)

def nmf(K: int, points: List[List[float]]) -> List[List[float]]:
    """Python wrapper function to calculate symnmf matrix by calling appropriate C module function.
    Args:
        K (int): number of clusters
        points (List[List[float]]): Datapoints used to calculate symnmf matrix
    Returns:
        List[List[float]]: Resultant symnmf matrix
    """
    norm_matrix = norm(points)
    H = initialize_H(norm_matrix, K)
    return symnmf_c.symnmf(H, norm_matrix)


def main():
    """Main function for SymNMF in Python, handles performing algorithm as per project specifications by calling
    appropriate wrapper functions based on user input, handles basic input verification.
    """
    args: argparse.Namespace = parse()
    K, goal, file_name = args.K, args.goal, args.file_name
    filepath = os.path.join(os.path.join(os.getcwd()), file_name)
    points = read_file(filepath)

    goals_mapping = {"symnmf": 0, "sym": 1, "ddg": 2, "norm": 3}
    if goal not in goals_mapping:
        print("An Error Has Occurred")
        return
    
    if goals_mapping[goal] == 0:
        try:
            K = float(K)
            if K != int(K):
                raise ValueError()
            K = int(K)
            
        except ValueError:
            print("An Error Has Occurred")
            return
        
        if K <= 1 or K >= len(points):
            print("An Error Has Occurred")
            return
        pretty_print(nmf(K, points))
    elif goals_mapping[goal] == 1:
        pretty_print(sym(points))
    elif goals_mapping[goal] == 2:
        pretty_print(diag(points))
    elif goals_mapping[goal] == 3:
        pretty_print(norm(points))



if __name__ == "__main__":
    main()

