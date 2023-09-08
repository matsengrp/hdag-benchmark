from sklearn.manifold import MDS, TSNE, Isomap
from sklearn.cluster import AgglomerativeClustering, SpectralClustering, DBSCAN
from sklearn.metrics import silhouette_score, silhouette_samples
import pickle
import historydag as hdag
import numpy as np
from collections import Counter

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

num_trees = 20  # Number of trees to downsample to (will grow e.g., 10 -> 50, 15 -> 80)

def main():
    # TODO:
    # - Load small DAG with >1000 trees
    # - Compute pairwise RF distance for each tree
    # - Use MDS on distance matrix
    #     - Could also try PCA or t-SNE
    # - Plot result for each tree
    #     - Could also highlight all the trees that are median. 

    print("Loading DAG...")
    dag = load_DAG()
    dag.summary()

    print("Loading pairwise distances...")
    # dists, is_med = get_pairwise_dist(dag)
    with open("dist_mat_complete.pkl", "rb") as f:
        dists, is_med = pickle.load(f)

    print("Clustering trees...")
    labels = cluster_label(dists, n_clusters=5)
    
    print("Computing sscores...")
    sscores = silhouette_score(dists, labels)
    sample_sscores = silhouette_samples(dists, labels)
    print("Silhoutte scores:", sscores)
    print("Median tree silhoutte score:", sample_sscores[is_med])

    plot_silhoutte_scores(dists, clustering_type="dbscan", embedding_type="tsne")
    # plot_silhoutte_scores(dists, is_med, linkage_method="single")
    # plot_silhoutte_scores(dists, is_med, linkage_method="average")
    
    # plot_cluster_embs(dists, is_med, n_clusters=5)
    # plot_cluster_embs(dists, is_med, n_clusters=4)
    # plot_cluster_embs(dists, is_med, n_clusters=6)
    

    # TODO:
    # - Examine properties of different tree clusters.
    #   - Node count, whether they group particular taxa together, etc.
    #   - Could the DAG be trimmed to contain a single cluster
    # - Silhoutte clustering analysis (compare to random groupings of trees?)
    # - Look at the agglomerative clustering dendrogram (i.e., what distance
    #   values between clusters are used as threshold).

def plot_silhoutte_scores(dists, clustering_type="agglomerative", embedding_type="tsne"):
    range_n_clusters = [5, 4, 6, 3]
    central_idxs, outlier_idxs = get_topologically_interesting_points(dists)

    for n_clusters in range_n_clusters:
        # Create a subplot with 1 row and 2 columns
        fig, (ax1, ax2) = plt.subplots(1, 2)
        fig.set_size_inches(10, 4.8)

        # The 1st subplot is the silhouette plot
        # The silhouette coefficient can range from -1, 1 but in this example all
        # lie within [-0.1, 1]
        ax1.set_xlim([-0.1, 1])
        
        # Initialize the clusterer with n_clusters value and a random generator
        # seed of 10 for reproducibility.
        if clustering_type == "agglomerative":
            clusterer = AgglomerativeClustering(metric='precomputed', linkage="complete", n_clusters=n_clusters)
            labels = clusterer.fit_predict(dists)
        elif clustering_type == "spectral":
            clusterer = SpectralClustering(affinity='precomputed', n_clusters=n_clusters)
            delta = 50
            labels = clusterer.fit_predict(np.exp(- dists ** 2 / (2. * delta ** 2)))
        elif clustering_type == "dbscan":
            clusterer = DBSCAN(metric='precomputed', eps=2, min_samples=20)
            labels = clusterer.fit_predict(dists)
            if -1 in labels:
                labels += 1
            n_clusters = len(np.unique(labels))

        # The (n_clusters+1)*10 is for inserting blank space between silhouette
        # plots of individual clusters, to demarcate them clearly.
        ax1.set_ylim([0, len(dists) + (n_clusters + 1) * 10])

        # The silhouette_score gives the average value for all the samples.
        # This gives a perspective into the density and separation of the formed
        # clusters
        silhouette_avg = silhouette_score(dists, labels)
        print(
            "For n_clusters =",
            n_clusters,
            "The average silhouette_score is :",
            silhouette_avg,
        )

        # Compute the silhouette scores for each sample
        sample_silhouette_values = silhouette_samples(dists, labels)

        y_lower = 10
        for i in range(n_clusters):
            # Aggregate the silhouette scores for samples belonging to
            # cluster i, and sort them
            ith_cluster_silhouette_values = sample_silhouette_values[labels == i]

            ith_cluster_silhouette_values.sort()

            size_cluster_i = ith_cluster_silhouette_values.shape[0]
            y_upper = y_lower + size_cluster_i

            color = cm.nipy_spectral(float(i) / n_clusters)
            ax1.fill_betweenx(
                np.arange(y_lower, y_upper),
                0,
                ith_cluster_silhouette_values,
                facecolor=color,
                edgecolor=color,
                alpha=0.7,
            )

            # Label the silhouette plots with their cluster numbers at the middle
            ax1.text(-0.05, y_lower + 0.5 * size_cluster_i, str(i))

            # Compute the new y_lower for next plot
            y_lower = y_upper + 10  # 10 for the 0 samples

        ax1.set_title("The silhouette plot for the various clusters.")
        ax1.set_xlabel("The silhouette coefficient values")
        ax1.set_ylabel("Cluster label")

        # The vertical line for average silhouette score of all the values
        ax1.axvline(x=silhouette_avg, color="red", linestyle="--")

        ax1.set_yticks([])  # Clear the yaxis labels / ticks
        ax1.set_xticks([-0.1, 0, 0.2, 0.4, 0.6, 0.8, 1])

        # 2nd Plot showing the actual clusters formed
        if embedding_type == "tsne":
            embedding = TSNE(n_components=2, metric="precomputed", init='random', random_state=1234)
        elif embedding_type == "isomap":
            embedding = Isomap(n_components=2, metric="precomputed", radius=4, n_neighbors=None)
        coords = embedding.fit_transform(dists)
        x = [coord[0] for coord in coords]
        y = [coord[1] for coord in coords]

        colors = cm.nipy_spectral(labels.astype(float) / n_clusters)
        ax2.scatter(
            x, y, marker=".", s=50, lw=0, alpha=0.7, c=colors, edgecolor="k"
        )
        ax2.scatter(x[central_idxs[0]], y[central_idxs[0]], c="red", s=75, marker="*", label="Central Trees")
        if len(central_idxs) > 1:
            for i in central_idxs[1:]:
                ax2.scatter(x[i], y[i], c="red", s=75, marker="*")

        ax2.scatter(x[outlier_idxs[0]], y[outlier_idxs[0]], c="blue", s=75, marker="*", label="Outlier Trees")
        if len(outlier_idxs) > 1:
            for i in outlier_idxs[1:]:
                ax2.scatter(x[i], y[i], c="blue", s=75, marker="*")
        ax2.legend()

        ax2.set_title("The visualization of the clustered data.")
        ax2.set_xlabel("Feature space for the 1st feature")
        ax2.set_ylabel("Feature space for the 2nd feature")


        plt.suptitle(
            "Silhouette analysis for Agglomerative clustering on sample data with n_clusters = %d"
            % n_clusters,
            fontsize=14,
            fontweight="bold",
        )

        plt.savefig(f"{clustering_type}_{embedding_type}_silhoutte_analysis_{n_clusters}.png")
        plt.clf()
        if clustering_type == "dbscan":
            return



def cluster_label(dist_matrix, n_clusters=5):
    clustering = AgglomerativeClustering(metric='precomputed', linkage="complete", n_clusters=n_clusters).fit(dist_matrix)
    return clustering.labels_


def load_DAG():
    with open('../historydag/sample_data/toy_trees.p', 'rb') as fh:
        ete_trees = pickle.load(fh)
    dag = hdag.history_dag_from_etes(ete_trees, ['sequence'])
    count = dag.count_histories()
    print(f"DAG contains {count} trees")

    # "Complete" the DAG, adding all allowable edges.
    dag.make_complete()
    dag.count_histories()  # 3431531

    # Show counts of trees with various parsimony scores.
    dag.hamming_parsimony_count()

    # "Trim" the DAG to make it only display minimum-weight trees.
    dag.trim_optimal_weight()
    # With default args, same as hamming_parsimony_count
    dag.weight_count()  # Counter({75: 45983})

    # "Collapse" the DAG, contracting zero-weight edges.
    dag.convert_to_collapsed()

    dag.weight_count()  # Counter({75: 1208})
    num_topologies = dag.count_topologies()  # 1054 unique topologies, ignoring internal labels
    print(f"DAG contains {num_topologies} topologies")

    # TODO: Try sampling x trees and building a dag from them...
    # dag.make_uniform()
    # samples = [dag.sample() for i in range(num_trees)]
    # dag = samples[0]
    # for other in samples[1:]:
    #     dag |= other

    # print(f"Downsampled DAG contains {dag.count_histories()} trees before trim")
    # dag.trim_optimal_weight()
    # print(f"Downsampled DAG contains {dag.count_histories()} trees after trim")

    return dag

def get_topologically_interesting_points(dists):
    u, idxs = np.unique(dists, axis=1, return_index=True)

    summed_rf = np.sum(dists[idxs], axis=1).reshape(-1)
    plt.hist(summed_rf, bins=25)
    plt.savefig("srf.png")
    plt.clf()

    max = np.max(summed_rf)
    min = np.min(summed_rf)
    print("Max SRF =", max)
    
    outlier_idxs = []
    for i, srf in enumerate(summed_rf):
        if srf >= max * (0.90):
            outlier_idxs.append(int(i))
    
    central_idxs = []
    for i, srf in enumerate(summed_rf):
        if srf <= min * (1.01):
            central_idxs.append(int(i))

    return np.array(central_idxs), np.array(outlier_idxs)

def get_pairwise_dist(dag):
    # Returns the pairwise distances and which trees are median trees in our dag
    n = dag.count_histories()
    dist = np.zeros((n, n))
    print("Creating tree list...")
    tree_list = [t for t in dag]
    for i, t1 in enumerate(tree_list):
        if i % 10 == 0:
            print(i)
        rf_dist = hdag.utils.make_rfdistance_countfuncs(t1)
        for j, t2 in enumerate(tree_list[i+1:], i+1):
            rf = t2.optimal_weight_annotate(**rf_dist)
            dist[i, j] = rf
            dist[j, i] = rf

    print("Finding medians")
    dag.trim_optimal_sum_rf_distance(dag)

    is_med = []
    num_medians = dag.count_histories()
    print(f"\tFound {num_medians}")
    for t in tree_list:
        new_dag = dag | t
        is_med.append(num_medians == new_dag.count_histories())

    assert sum(is_med) == num_medians

    with open("dist_mat.pkl", "wb") as f:
        pickle.dump((dist, is_med), f)

    return dist, is_med

def plot_embs(dist_matrix, is_med):
    # embedding = MDS(n_components=2, dissimilarity='precomputed')
    embedding = TSNE(n_components=2, metric="precomputed", init='random')
    coords = embedding.fit_transform(dist_matrix)
    x = [coord[0] for coord in coords]
    y = [coord[1] for coord in coords]

    med_idxs = [i for i, med in enumerate(is_med) if med]   # I think these are all the same
    med_idx = med_idxs[0]

    plt.scatter(x, y, alpha=0.5, c=dist_matrix[med_idx])
    plt.colorbar()
    plt.scatter(x[med_idx], y[med_idx], c="red", marker="*", label="Median Tree")
    plt.title("t-SNE Projection of Pairwise RF Distances for 1k MP trees")
    plt.legend()
    plt.savefig("tree_embs.pdf")

def plot_cluster_embs(dist_matrix, is_med, n_clusters=5):
    # embedding = MDS(n_components=2, dissimilarity='precomputed')
    embedding = TSNE(n_components=2, metric="precomputed", init='random')
    coords = embedding.fit_transform(dist_matrix)
    x = [coord[0] for coord in coords]
    y = [coord[1] for coord in coords]

    med_idxs = [i for i, med in enumerate(is_med) if med]   # I think these are all the same
    med_idx = med_idxs[0]

    clustering = AgglomerativeClustering(metric='precomputed', linkage="complete", n_clusters=n_clusters).fit(dist_matrix)
    # clustering = AgglomerativeClustering(metric='precomputed', linkage="complete", distance_threshold=15, n_clusters=None).fit(dist_matrix)
    hues = clustering.labels_

    cluster_count = Counter()
    for hue in hues:
        cluster_count[hue] += 1

    print(cluster_count)

    plt.scatter(x, y, alpha=0.5, c=hues)
    plt.scatter(x[med_idx], y[med_idx], c="red", marker="*", label="Median Tree")
    # plt.colorbar()
    plt.title("t-SNE Projection of Pairwise RF Distances for 1k MP trees")
    plt.legend()
    plt.savefig(f"tree_embs_clusters_{n_clusters}.pdf")
    plt.clf()


def plot_med_dist_hist(dists, is_med):

    med_idxs = [i for i, med in enumerate(is_med) if med]   # I think these are all the same
    med_idx = med_idxs[0]

    num_bins = int(np.max(dists[med_idx]) - np.min(dists[med_idx]))
    print(num_bins)
    plt.hist(dists[med_idx], bins=num_bins)
    plt.ylabel("Count")
    plt.xlabel("RF Distance to Median Tree")
    plt.title("Distribution of Pairwise RF Distances")
    plt.savefig("tree_med_dist_hist.pdf")





if __name__ == "__main__":
    main()