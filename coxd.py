__author__ = 'maximkuleshov'

from itertools import combinations
from statistics import median, mean


def get_coex_dist(genes, engine, id_gene):
    coex_median = 0
    dist_median = 0
    if len(genes) > 1:
        pairs = combinations(genes, 2)
        wrap = 'SELECT * FROM coex_dist NATURAL JOIN (SELECT {0}) t;'
        template = ''' {0} AS gene1, {1} AS gene2 '''

        query_pairs = []
        for pair in pairs:
            g1, g2 = sorted(pair)
            if (g1 in id_gene) and (g2 in id_gene):
                query_pairs.append(template.format(id_gene[g1], id_gene[g2]))

        res = []
        if query_pairs:
            query = wrap.format(' UNION ALL SELECT '.join(query_pairs))
            res = list(engine.execute(query))

        try:
            coex_median = median(line[3] for line in res)
            dist_median = mean(line[4] for line in res)
        except:
            return '{0}\t{1}'.format(coex_median, dist_median)

    return '{0}\t{1}'.format(coex_median, dist_median)


def main():
    return None


if __name__ == '__main__':
    main()
