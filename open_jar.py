__author__ = 'maximkuleshov'

import dill as pickle
import matplotlib.pyplot as plt
from matplotlib import rc, rcParams, rcParamsDefault
from math import log10


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


def plot_ranks(levels, length):
    plt.axes([0.15, 0.68, 0.8, 0.27])
    plt.xlim([0, length])
    plt.plot(levels)
    plt.ylabel('rank level, %', fontsize=10, labelpad=15)


def plot_corr(indices, length):
    x = sorted(indices.keys())
    y = [indices[key][0] for key in x]

    plt.axes([0.15, 0.37, 0.8, 0.27])
    plt.xlim([0, length])
    plt.scatter(x, y)
    plt.ylabel('-log10(p-value)*corr', fontsize=10, labelpad=15)


def plot_dist(indices, length):
    x = sorted(indices.keys())
    y = [indices[key][1] for key in x]

    plt.axes([0.15, 0.06, 0.8, 0.27])
    plt.xlim([0, length])
    plt.scatter(x, y)
    plt.ylabel('-log10(p-value)*(3-dist)', fontsize=10, labelpad=15)


def main():
    chea = 'ChEA_2016'
    direction = 'up'
    myc0 = pickle.load(open('%s_%s_pval.05.pickle' % (chea, direction), 'rb'))[1]['MYC']

    plot_data = res2plot(myc0['gene:1007'])

    plt.figure(figsize=(6, 8))
    plot_ranks(plot_data[0], plot_data[2])
    plot_corr(plot_data[1], plot_data[2])
    plot_dist(plot_data[1], plot_data[2])

    plt.xlabel('position')
    plt.show()

    return



if __name__ == '__main__':
    main()