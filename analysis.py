from symnmf import symnmf as snmf
from kmeans import run_kmeans as kmeans, load_points_from_file as load_data, euclidian_distance as dist
from sklearn.metrics import silhouette_score
import sys

if __name__ == "__main__":

    num_clusters, data_file, convergence_threshold, max_iterations = 0, None, 0.0001, 300
    num_clusters, data_file = int(sys.argv[1]), sys.argv[2]

    data_matrix = load_data(data_file)

    factorization_matrix = snmf('symnmf', num_clusters, data_file)
    
    
    initial_centroids = kmeans(data_matrix, num_clusters, max_iterations, convergence_threshold)

    # Assign clusters for the factorization method (SymNMF)
    symnmf_clusters = [(max([(cl) for cl in range(len(factorization_matrix[0]))], key=lambda x: factorization_matrix[i][x]))
                       for i in range(len(factorization_matrix))]

    # Assign clusters for the k-means method
    kmeans_clusters = [(min([k for k in range(len(initial_centroids))], key=lambda k: dist(
        data_matrix[i], initial_centroids[k]))) for i in range(len(data_matrix))]

    print("SymNMF Silhouette:", '%.4f' % silhouette_score(data_matrix, symnmf_clusters))
    print("K-Means Silhouette:", '%.4f' % silhouette_score(data_matrix, kmeans_clusters))