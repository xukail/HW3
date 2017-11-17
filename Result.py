from LookForORF import look_for_orfs, Viterbi
# from sklearn.metrics import roc_curve, auc
import collections
import math
import matplotlib.pyplot as plt

def main():
    print '(a) The language I used is Python 2.7'
    print '(b) Basic info about computer: (Run on VMVare Intel Core i7-6700HQ CPU 2.60GHz * 2 with 3 GiB of RAM)'
    print 'the total CPU time taken by your algorithm for the first 9 Viterbi iterations: '

    res, genome = look_for_orfs()
    # print 'finished initialiing'
    Viterbi(genome, 10)

if __name__ == "__main__":
    main()
    print 'Done.'