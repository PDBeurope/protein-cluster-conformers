"""
Performs agglomerative clustering on the sum-based scores generated from CA distance
difference matrices.

Script to plot a dendrogram of the hierarchical clustering results. Can either load the
results from a pre-saved CSV of the clustering output from cluster.py or work on the
benchmark dataset.

Functions included herein to generate a swarm plots from a given dir of pre-generated
distance difference maps.
"""

# Third-party modules
import seaborn as sns
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from numpy import column_stack, ndarray, zeros
from scipy.cluster.hierarchy import dendrogram
from sklearn.cluster import AgglomerativeClustering

from cluster_conformers.utils.linear_algebra_utils import upper_triangle


def make_linkage_matx(model: AgglomerativeClustering) -> "ndarray[any, float]":
    """
    Linkage matrix returned from a SKLearn model for agglomerative clustering.
    Returned: Numpy ndarray

    :param model: SKLearn agglomerative clustering model, fitted to some score matrix.
    :type model: sklearn.cluster.AgglomerativeClustering
    :return: Matrix of descriptors for the clustering results from the SKLearn model.
    :rtype: np.ndarray[any, float]
    """

    counts = zeros(model.children_.shape[0])
    n_samples = len(model.labels_)
    for i, merge in enumerate(model.children_):
        current_count = 0
        for child_idx in merge:
            if child_idx < n_samples:
                current_count += 1  # leaf node
            else:
                current_count += counts[child_idx - n_samples]
        counts[i] = current_count

    linkage_matrix = column_stack([model.children_, model.distances_, counts]).astype(
        float
    )

    return linkage_matrix


def cluster_agglomerative(
    score_matx: "ndarray[any, float]", cutoff: float = None
) -> AgglomerativeClustering:
    """
    Performs the agglormerative (bottom up) clustering algorithm on a given score
    matrix.

    :param score_matx: N*N-square, symmetric matrix of scores used for clustering.
    :type score_matx: ndarray[any, float]
    :param cutoff: Threshold below which to begin defining groups as clusters, defaults
        to None
    :type cutoff: float, optional
    :return: Scikit-learn Agglomerative model fitted to input data.
    :rtype: sklearn.cluster.AgglomerativeClustering
    """

    # Get the cutoff
    model = AgglomerativeClustering(
        n_clusters=None,  # Force function to fix number of clusters
        distance_threshold=0,
        metric="precomputed",  # Force to use score matrix
        linkage="average",  # UPGMA
        compute_distances=True,
    )

    # Repeat but cluster based on parsed cutoff
    if cutoff:

        # Fit score matrix to model
        model = model.fit(score_matx)
        linkage_matx = make_linkage_matx(model)

        cutoff *= linkage_matx[:, 2].max()

        model = AgglomerativeClustering(
            n_clusters=None,
            distance_threshold=cutoff,
            metric="precomputed",
            linkage="average",
            compute_distances=True,
        )

    # Fit score matrix to model
    model = model.fit(score_matx)

    return model


def plot_dendrogram(
    unp: str, axis, linkage_matrix: ndarray = None, cutoff: float = None, **kwargs
) -> "tuple[Figure, Axes]":
    """
    Create linkage matrix from SKLearn model and plot the dendrogram of nodes. Applies
    this to parsed Matplotlib axis, inplace.

    :param unp: UniProt accession
    :type unp: str
    :param linkage_matrix: Matrix of clustering information (labels and distances)
        derived from the fitted agglomerative clustering model.
    :type linkage_matrix: np.ndarray
    :param cutoff: Distance threshold below which clusters were defined, has no affect
        on cluster results and is simply used to plot horizontal line for clarity.
        Should be 0-1, defaults to None.
    :param kwargs: Same parameters as accepted by scipy.cluster.hierarchy.dendrogram()
    :type kwargs: Any
    :type cutoff: float, optional
    :return: Matplotlib figure and axis.
    :rtype: tuple(matplotlib.figure.Figure, matplotlib.axes.Axes)
    """

    # Plot the corresponding dendrogram
    dendrogram_plot = dendrogram(linkage_matrix, ax=axis, **kwargs)

    hline_props = {
        "colors": ["black"],
        "linestyles": ["dashed"],
        "linewidths": 1,
        "alpha": 0.5,
    }

    # Add horizontal line where cutoff is placed
    if cutoff:
        max_parent = linkage_matrix[:, 2].max()
        axis_xlimits = axis.get_xlim()
        axis.hlines(
            y=max_parent * cutoff,
            xmin=axis_xlimits[0],
            xmax=axis_xlimits[1],
            **hline_props,
        )

    axis.set_title(f"Agglomerative clustering results: {unp}", fontweight="bold")
    axis.set_ylabel("Score (\u212B)")

    # Set the y-axis limits
    rock_bottom = -max_parent * 0.025
    axis.set_ylim(rock_bottom, max_parent * 1.05)

    # Add vlines for each leaf below the x-axis
    axis.vlines(
        x=axis.get_xticks(),
        ymin=rock_bottom,
        ymax=0,
        colors=dendrogram_plot["leaves_color_list"],
        linestyles=["solid"],
        linewidths=1,
        # alpha=0.5,
    )

    axis.hlines(y=0, xmin=axis_xlimits[0], xmax=axis_xlimits[1], **hline_props)

    del dendrogram_plot


def plot_swarmplot(scores: ndarray, axis) -> "tuple[Figure, Axes]":
    """Creates a swarm plot of non-overlapping data points for a given list of data. The
    values of the data will correspond to their y-values. Their position along the
    x-axis is irrelevant as they're all identical.

    :param scores: Matrix of GLOCON scores.
    :type scores: Iterable
    :param axis: Matplotlib axis to plot the swarm plot on.
    :type axis: Axes
    :return: Figure and axis objects containing the plotted swarm plot
    :rtype: tuple(matplotlib.figure.Figure, matplotlib.axes.Axes)
    """
    # Plot the data
    sns.swarmplot(data=upper_triangle(scores), ax=axis, size=5)  # vmax=max_dist ,

    # Add some formatting
    axis.set_ylabel("GLOCON score (\u212B)")
