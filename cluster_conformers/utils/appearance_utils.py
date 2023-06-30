"""
Functions for modifying data visualisation behaviour and calculation progression.
"""

from matplotlib import pyplot as plt


def init_plot_appearance():
    """
    Setup of default appearance for any plot
    """

    # Figure parameters and formatting
    plt.style.use("seaborn-v0_8-colorblind")
    plt.style.use("fast")
