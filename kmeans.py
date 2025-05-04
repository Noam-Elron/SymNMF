from multiprocessing import Value
import os
from typing import List
import math
import argparse
import sys

epsilon = 0.0001

class Point:
    """Class Wrapper for a Datapoint
    """
    def __init__(self, point=None):
        if point != None:
            self.vals = point
        else:
            self.vals = []
        
    def push(self, val: float):
        """Standard push operation to add coordinate to datapoint

        Args:
            val (float): coordinate to be added to datapoint
        """
        self.vals.append(val)

    def euclidean_distance(self, other: "Point") -> float:
        """Calculates euclidean distance between current point and another

        Args:
            other (Point): Datapoint we are calculating euclidean distance from

        Returns:
            float: Euclidean distance between the two points
        """
        total = 0
        for i in range(len(self.vals)):
            total += (self.vals[i] - other.vals[i]) * (self.vals[i] - other.vals[i])
        return math.sqrt(total)

    def __add__(self, other: "Point") -> "Point":
        """Adds current and other point together, coordinate by coordinate and returns the result. (Does not mutate existing points) 

        Args:
            other (Point): Datapoint we're adding to our own

        Returns:
            Point: Resultant datapoint from add operation
        """
        new_point = Point()
        for i in range(len(self.vals)):
            x_i = self.vals[i] + other.vals[i]
            new_point.push(x_i)
        return new_point

    def __truediv__(self, num: int) -> "Point":
        """Divides current point by an integer, coordinate by coordinate and returns the result. (Does not mutate existing point) 

        Args:
            num (int): Integer to divide current point by

        Returns:
            Point: _Resultant datapoint from division operation
        """
        assert type(num) == int

        new_point = Point()
        for i in range(len(self.vals)):
            x_i = self.vals[i] / num
            new_point.push(x_i)
        return new_point
    
    def __repr__(self):
        return f"Point({self.vals})"
    
    def __str__(self) -> str:
        """String formatter to print datapoints as per course specifications

        Returns:
            str: String representation of current point
        """
        formatted_floats = list(map(lambda x: f'{x:.4f}', self.vals))
        result = ', '.join(formatted_floats)

        return result


class Cluster:
    """Class wrapper for Centroids/Cluster
    """
    def __init__(self, centroid):
        self.centroid: Point = centroid        
        self.total_points: Point = Point([0 for i in range(len(self.centroid.vals))]) 
        self.num_points: int = 0

    def update_centroid(self, point: Point):
        """Modifies current cluster/centroid with an additional datapoint

        Args:
            point (Point): Datapoint we're adding to cluster
        """
        self.total_points += point
        self.num_points += 1

    def finalize_next_centroid_pos(self) -> float:
        """Perform cluster operation to find new current cluster centroid. Returns distance of previous centroid to new one

        Returns:
            float: distance of previous centroid to new one
        """
        if self.num_points == 0:
            return 0
        prev_centroid = self.centroid
        self.centroid = self.total_points / self.num_points
        self.num_points = 0
        self.total_points = Point([0 for i in range(len(self.centroid.vals))])
        return prev_centroid.euclidean_distance(self.centroid)
    
    def __str__(self):
        return str(self.centroid)


class ClusterManager:
    """Wrapper for class that wraps helper functions for Cluster/Centroid management
    """
    def __init__(self):
        self.clusters: List[Cluster] = []

    def closest_cluster(self, point: Point) -> int:
        """Finds and returns ID (index) of closest centroid/cluster

        Args:
            point (Point): Point we want to find closest cluster to

        Returns:
            int: Index of closest cluster to point
        """
        closest_id = 0
        closest_dist = point.euclidean_distance(self.clusters[0].centroid)
        for i in range(1, len(self.clusters)):
            dist = point.euclidean_distance(self.clusters[i].centroid)
            if dist < closest_dist:
                closest_id = i
                closest_dist = dist
        return closest_id

    def update_next_centroid_position(self, centroid_id: int, point: Point):
        """Helper/Wrapper function to update cluster's centroid by adding a new point to it

        Args:
            centroid_id (int): ID of centroid we want to update
            point (Point): Point we want to add to cluster
        """
        self.clusters[centroid_id].update_centroid(point)

    def push(self, cluster: Cluster) -> None:
        """Adds new cluster to ClusterManager

        Args:
            cluster (Cluster): New cluster we're adding
        """
        self.clusters.append(cluster)

    def __iter__(self):
        return iter(self.clusters)
    
    def clusters_to_list(self):
        """Helper function to be able to easily pass info at end of KMeans operation to standard python output for easy usage

        Returns:
            List[List[float]]: 2D Array of clusters
        """
        clusters = []
        for cluster in self:
            clusters.append(cluster.centroid.vals)
        return clusters

def read_points(input_data: str) -> List[Point]:
    """Function to read datapoints from file
    Args:
        input_data (str): filepath to .txt file containing valid datapoints
    Returns:
        List[Point]: 2D Array, each element is a point
    """
    points: List[Point] = []
    with open(input_data, "r", encoding="utf-8") as input:
        for line in input:
            line = line.strip()
            line = line.split(",")
            line = list(map(float, line))
            point = Point(line)
            points.append(point)
    return points

def kmeans(K: int, iter: int, input_data: str) -> List[List[float]]:
    """KMeans Wrapper for entire functionality

    Args:
        K (int): Number of clusters
        iter (int): Maximum number of iterations in KMeans
        input_data (str): filepath to .txt file containing the datapoints

    Returns:
        List[List[float]]: 2D array of clusters resulting from KMeans algorithm
    """
    if iter >= 1000 or iter <= 1:
        print("Invalid maximum iteration!")
        exit(1)
    
    points: List[Point] = read_points(input_data)
            
    if K <= 1 or K >= len(points):
        print("Invalid number of clusters!")
        exit(1)
    
    manager = ClusterManager()

    for i in range(K):
        cluster = Cluster(points[i])
        manager.push(cluster)

    for i in range(iter):
        for point in points:
            nearest_cluster = manager.closest_cluster(point)
            manager.update_next_centroid_position(nearest_cluster, point)
        
        converge = True
        for cluster in manager: 
            delta = cluster.finalize_next_centroid_pos()
            if delta > epsilon:
                converge = False
        
        if converge:
            break
    
    return manager.clusters_to_list()

def parse() -> argparse.Namespace:
    """Parsing Function for user input
    Returns:
        argparse.Namespace: Parsing object to pass user input to main program
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('K', type=str)
    parser.add_argument('iter', type=str, nargs='?', default="200")
    parser.add_argument('input_file', type=str)

    if len(sys.argv) < 3 or len(sys.argv) > 4:
        print("An Error Has Occurred")
        exit(1)
        
    return parser.parse_args()

if __name__ == "__main__":
    args: argparse.Namespace = parse()
    K, iterations, filepath = args.K, args.iter, args.input_file
    
    try:
        K = float(K)
        if K != int(K):
            raise ValueError()
        K = int(K) 
    except ValueError:
        print("Invalid number of clusters!")
        exit(1)

    try:
        iterations = float(iterations)
        if iterations != int(iterations):
            raise ValueError()
        iterations = int(iterations)
    except ValueError:
        print("Invalid maximum iteration!")
        exit(1)

    try:
        filepath = os.path.join(os.getcwd(), filepath)
        kmeans_output = kmeans(K, iterations, filepath)
        for cluster_output in kmeans_output:
            formatted_floats = list(map(lambda x: f'{x:.4f}', cluster_output))
            result = ','.join(formatted_floats)
            print(result)
    except ValueError:
        print("An Error Has Occurred")


