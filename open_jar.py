__author__ = 'maximkuleshov'

import dill as pickle

def main():
    chea = 'ChEA_2016'
    direction = 'up'
    jar = pickle.load(open('%s_%s_pval.05.pickle' % (chea, direction), 'rb'))
    return None


if __name__ == '__main__':
    main()