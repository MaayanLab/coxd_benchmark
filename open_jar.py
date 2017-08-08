__author__ = 'maximkuleshov'

import dill as pickle
import matplotlib.pyplot as plt
from matplotlib import rc, rcParams, rcParamsDefault
from math import log10, sqrt


def res2plot(result):
    inc = 100/len(result)
    length = result[-1][0]+10
    current_level = 0
    rank_levels = []
    indices = {i[0]-1: [-log10(i[2])*i[9], -log10(i[2]*(3-i[10]))] for i in result}

    for pos in range(length):
        if pos+1 in indices:
            current_level += inc
        rank_levels.append(current_level)
    return [rank_levels, indices, length]


def plot_ranks(levels, length, ax):
    ax.set_xlim([0, length])
    ax.plot(levels)
    # ax.set_ylabel('rank level, %')


def plot_corr(indices, length, ax):
    x = sorted(indices.keys())
    y = [indices[key][0] for key in x]

    ax.set_xlim([0, length])
    ax.scatter(x, y)
    # ax.set_ylabel('-log10(p-value)*corr')


def plot_dist(indices, length, ax):
    x = sorted(indices.keys())
    y = [indices[key][1] for key in x]

    ax.set_xlim([0, length])
    ax.scatter(x, y)
    # ax.set_ylabel('-log10(p-value)*(3-dist)')
    ax.set_xlabel('position')


def main():
    chea = 'ChEA_2016'
    direction = 'up'
    myc = pickle.load(open('%s_%s_pval.05.pickle' % (chea, direction), 'rb'))[1]['MYC']

    fig_width = 12
    fig_height = 8
    rows = int(sqrt(len(myc.keys())))
    cols = rows+1

    fig, axes = plt.subplots(nrows=rows, ncols=cols, figsize=(fig_width, fig_height))

    for i, key in enumerate(sorted(myc.keys())):
        plot_data = res2plot(myc[key])
        plot_ranks(plot_data[0], plot_data[2], axes.flat[i])
        plot_corr(plot_data[1], plot_data[2], axes.flat[i])
        plot_dist(plot_data[1], plot_data[2], axes.flat[i])

    plt.tight_layout()
    plt.show()



    return



if __name__ == '__main__':
    main()