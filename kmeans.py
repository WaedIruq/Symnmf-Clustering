import math
import sys


def euclidian_distance(x, y):
    """Compute the Euclidean distance between two points."""
    return math.sqrt(sum((x_i - y_i) ** 2 for x_i, y_i in zip(x, y)))


def group_points_by_closest_center(data_points, centers):
    """Assign each point to the nearest center."""
    clusters = [[] for _ in range(len(centers))]
    for point in data_points:
        closest_center_idx = min(range(len(centers)),
                                  key=lambda idx: euclidian_distance(point, centers[idx]))
        clusters[closest_center_idx].append(point)
    return clusters


def compute_new_centers(groups):
    """Calculate the new centers by averaging each group."""
    new_centers = []
    for group in groups:
        new_center = [sum(dim_value) / len(group) for dim_value in zip(*group)]
        new_centers.append(new_center)
    return new_centers


def run_kmeans(data_points, num_clusters, max_iterations, convergence_criteria):
    """Execute the K-means clustering algorithm."""
    centers = data_points[:num_clusters]
    for _ in range(max_iterations):
        previous_centers = centers[:]
        clusters = group_points_by_closest_center(data_points, centers)
        centers = compute_new_centers(clusters)
        
        # Check if the centers have converged
        if max(euclidian_distance(previous_centers[i], centers[i]) for i in range(num_clusters)) < convergence_criteria:
            break
    return centers


def display_clusters_centers(centers):
    """Print the final centers of the clusters."""
    for center in centers:
        print(','.join([f'{dim:.4f}' for dim in center]))
    print()


def load_points_from_file(file_path):
    """Read the data points from the file."""
    points = []
    with open(file_path, 'r') as f:
        for line in f:
            points.append([float(value) for value in line.strip().split(',')])
    return points


if __name__ == "__main__":
    arguments = sys.argv
    clusters_count = int(arguments[1])
    input_file_path = arguments[-1]
    max_iter = 200 if len(arguments) == 3 else int(arguments[2])
    epsilon_threshold = 0.001
    data_points = load_points_from_file(input_file_path)

    if not (1 < clusters_count < len(data_points)):
        print("Invalid number of clusters!")
    if not (1 < max_iter < 1000):
        print("Invalid maximum iterations!")

    if 1 < clusters_count < len(data_points) and 1 < max_iter < 1000:
        display_clusters_centers(run_kmeans(data_points, clusters_count, max_iter, epsilon_threshold))