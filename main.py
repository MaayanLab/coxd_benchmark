__author__ = 'maximkuleshov'

import json
import os.path
import pickle
from collections import defaultdict
from operator import itemgetter
from time import sleep
from joblib import delayed, Parallel

import requests
from retrying import retry
from sqlalchemy import create_engine

from coxd import get_coex_dist

ENRICHR_URL = 'http://amp.pharm.mssm.edu/Enrichr/'


@retry
def get_enrichr_results(gene_set_library, genelist, description):
    """ Sends a gene list to Enrichr, gets result's id and retrieves results in JSON format  """
    addlist_url = ENRICHR_URL + 'addList'
    payload = {
        'list': (None, genelist),
        'description': (None, description)
    }

    response = requests.post(addlist_url, files=payload)
    if not response.ok:
        raise Exception('Error analyzing gene list')
    sleep(1)
    data = json.loads(response.text)

    enrich_url = ENRICHR_URL + '/enrich'
    query_string = '?userListId=%s&backgroundType=%s'
    user_list_id = data['userListId']
    response = requests.get(enrich_url + query_string % (user_list_id, gene_set_library))
    sleep(1)
    return json.loads(response.text)


def parse_gmt(gmt, direction):
    """ Filter CREEDS by direction and collapse terms for similar TFs into dictionary """
    tfs = defaultdict(list)
    for line in gmt:
        term, desc, *genes = line.strip().split('\t')
        genes = [gene.split(',')[0] for gene in genes]
        if term.split(sep='-')[1] == direction:
            tf = term.split(sep='-')[0].upper()
            tfs[tf].append(genes)
    return tfs


def filter_library(ref, lib):
    """ Calculates intersection between CREEDS terms and ChEA terms and returns corresponding records from CREEDS """
    lib_keys = set(line.split()[0].split(sep='_')[0] for line in lib)
    filtered_keys = set(ref.keys()).intersection(lib_keys)
    filtered_ref = []
    for key in sorted(filtered_keys):
        filtered_ref.extend([[key, ref[key][pos]] for pos in range(len(ref[key]))])
    return sorted(filtered_ref)


def map_tf(tf, res, ref):
    """ Maps TF's position(s) to the ChEA libraey sized histogram """
    indices = [i for i, x in enumerate(res) if x == tf]
    for index in indices:
        ref[index] += 1
    return ref


def main():
    libraries = [['ChEA_2016', 645]]
    dirs = ['up', 'dn']
    gmt_file = 'single_gene_perturbations-v1.0.gmt'
    engine = create_engine(
        'mysql+pymysql://coexpression_and_distance:systemsbiology@amp.pharm.mssm.edu/coexpression_and_distance')
    connection = engine.connect()

    # Create id to gene mapping table and keep it in memory
    query = list(engine.execute("SELECT * FROM genes;"))
    ids = [gene[0] for gene in query]
    genes = [gene[1] for gene in query]
    id_gene = dict(zip(genes, ids))

    for direction in dirs:
        for lib_data in libraries:
            chea, chea_size = lib_data
            creeds = parse_gmt(open(gmt_file, 'r').readlines(), direction)
            lib_file = open('%s.gmt' % chea, 'r').readlines()
            creeds = filter_library(creeds, lib_file)

            if os.path.isfile('%s_%s_pval.05.pickle' % (chea, direction)):
                jar = pickle.load(open('%s_%s_pval.05.pickle' % (chea, direction), 'rb'))
                start_pos, pval_hist, adj_pval_hist = jar
            else:
                start_pos = 0
                # Create a ChEA library sized histogram
                pval_hist = [0] * chea_size
                adj_pval_hist = [0] * chea_size

            for pos, line in enumerate(creeds[start_pos:]):
                print('{0}/{1}'.format(start_pos + pos, len(creeds)))
                key, genes = line
                if not genes:
                    continue
                results = []

                # Filter results by p-value and tf
                data = sorted([i for i in get_enrichr_results(chea, '\n'.join(genes), '')[chea] if (i[2] <= 0.05)], key=itemgetter(2))
                data = [j for j in [[i] + line[1:] for i, line in enumerate(data)] if j[1].split('_')[0].upper() == key]

                if data:
                    for res in data:
                        overlap = res[5]
                        coex_dist = [float(i) for i in get_coex_dist(overlap, engine, id_gene).split('\t')]
                        res.extend(coex_dist)
                        results.append(res)

                    # Sort tfs from result
                    # 0 - Index, 1 - Term
                    # 2 - P-value
                    # 3 - Z-score
                    # 4 - Combined score
                    # 5 - Genes
                    # 6 - Adjusted p-value
                    # 7 - Old p-value
                    # 8 - Old adjusted p-value
                    # 9 - Coexpression
                    # 10 - Distance
                    s_pval = [line[0] for line in sorted(results, key=itemgetter(2))]
                    pval_hist = map_tf(key, s_pval, pval_hist)

                    s_adj_pval = [line[0] for line in sorted(results, key=itemgetter(6)) if line[6] <= 0.05]
                    adj_pval_hist = map_tf(key, s_adj_pval, adj_pval_hist)

                    status = [start_pos + pos + 1, pval_hist, adj_pval_hist]
                    pickle.dump(status, open('%s_%s_pval.05.pickle' % (chea, direction), 'wb'))
        connection.close()
        engine.dispose()
    return None


if __name__ == '__main__':
    main()
