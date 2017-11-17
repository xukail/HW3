import math
import time

def look_for_orfs():
    valid_letters = ['A', 'C', 'T', 'G']
    genome_FASTA_filename = 'GCF_000091665.1_ASM9166v1_genomic'
    f = open('data/' + genome_FASTA_filename + '.fna', 'r')
    f.readline()
    genome_list = list(f.read().upper())
    f.close()

    ommit = None
    i = 0
    while i < len(genome_list):
        if genome_list[i] == '>':
            ommit = i
            break
        if genome_list[i] == '\n':
            genome_list.pop(i)
            i -= 1
        elif genome_list[i] not in valid_letters:
            genome_list[i] = 'T'
        i += 1

    genome = ''.join(genome_list[0:ommit])

    res = None
    # res = look_for_stops(genome)
    # print res
    return res, genome

def look_for_stops(genome):
    stop_codons = ['TAA', 'TAG', 'TGA']
    stops_lists = []
    len_of_genome = len(genome)
    for x in range(0, 3):
        i = x;
        stops = []
        while i < len_of_genome:
            triplet = genome[i : i + 3]
            if triplet in stop_codons:
                stops.append(i + 3)
            i += 3;
        stops_lists.append(stops)
    return stops_lists

def Viterbi(genome, num_iter):
    """
        Run Viterbi algorithm on a sequence of genome for num_iter times with training
    """
    e1 = {'A':0.25, 'C':0.25, 'G':0.25, 'T':0.25}
    e2 = {'A':0.20, 'C':0.30, 'G':0.30, 'T':0.20}
    a1 = {1:0.9999, 2:0.0001}
    a2 = {1:0.01, 2:0.99}
    e = {1:e1, 2:e2}
    a = {1:a1, 2:a2}
    begin1 = 0.9999
    begin2 = 0.0001

    start_time = time.time()

    for n in range(num_iter):
        print 'Iteration ', (n + 1)
        v_k_i_dict = {1:math.log(begin1 * e[1][genome[0:1]]), 2:math.log(begin2 * e[2][genome[0:1]])}
        v_list = [v_k_i_dict.copy()]

        states = []

        k = 1
        for i in range(1, len(genome)):
            x_i_plus_1 = genome[i:i + 1]
            v_l_i_plus_1_dict = {}
            for l in range(1, 3):
                if v_k_i_dict[1] + math.log(a[1][l]) < v_k_i_dict[2] + math.log(a[2][l]):
                    k = 2
                else:
                    k = 1
                v_k_i = v_k_i_dict[k]

                v_l_i_plus_1_dict[l] = v(e, a, l, x_i_plus_1, k, v_k_i)
            v_k_i_dict = v_l_i_plus_1_dict.copy()
            v_list.append(v_k_i_dict.copy())

        num_hits = 5
        if n == num_iter - 1:
            num_hits = float('infinity')
        # print 'Start to traceback'
        begin1, log_prob, hits = traceback_Viterbi(states, v_list, a, e, genome, num_hits)
        begin2 = 1.0 - begin1
        print '\t(a) The HMM emission/transition parameters: '
        print '\t Emission: '
        print '\t\t State 1: ', e[1]
        print '\t\t State 2: ', e[2]
        print '\t Transitions\t: State1 \tState2'
        print '\t Begin: \t', begin1, ' \t', begin2
        print '\t State1: \t', a[1]
        print '\t State2: \t', a[2]
        print ''
        print '\t(b) the log probability (natural log, base-e) of the overall Viterbi path.'
        print '\t\t', log_prob
        print ''
        print '\t(c) the total number of \"hits\" found, where a hit is (contiguous) subsequence assigned to state 2 in the Viterbi path:'
        print '\t\t', len(hits) / 2
        print '\t(d) the lengths and locations (starting and ending positions) of the first k \"hits.\"'
        i = 0
        res_hits = []
        while i < len(hits) and num_hits > i / 2:
            res_hits.append(tuple([hits[i], hits[i + 1]]))
            i += 2
        print '\t\t', res_hits
            
    # print("--- %s seconds ---" % (time.time() - start_time))



def v(e, a, l, x_i_plus_1, k, v_k_i):
    """
        The Viterbi algorithm calculate the log(based e) probability of a sequence ending in
        l state and x_i_plus_1 letter.
    """
    return math.log(e[l][x_i_plus_1]) + v_k_i + math.log(a[k][l])


def traceback_Viterbi(states, v_list, a, e, genome, num_hits):
    """
        Traceback to predict hidden states and calculate new emmisions and transitions
    """
    e1_count = {'A':0, 'C':0, 'G':0, 'T':0}
    e2_count = {'A':0, 'C':0, 'G':0, 'T':0}
    a1_count = {1:0, 2:0}
    a2_count = {1:0, 2:0}
    e_count = {1: e1_count, 2: e2_count}

    i = len(v_list) - 1
    l = 1
    if v_list[i][2] > v_list[i][1]:
        l = 2

    v_l_i_plus_1 = v_list[i][l]
    e_count[l][genome[i:i+1]] += 1
    # states.append(l)

    state_2 = []  # list of state2 starts and ends
    
    while i > 0:
        v_k_i_dict = v_list[i - 1]
        if v_k_i_dict[1] + math.log(a[1][l]) > v_k_i_dict[2] + math.log(a[2][l]):
            if l == 2:
                state_2.append(i + 1)
                a1_count[2] += 1
            else:
                a1_count[1] += 1
            l = 1
        else:
            if l == 1:
                state_2.append(i)
                a2_count[1] += 1
            else:
                a2_count[2] += 1
            l = 2
        i -= 1
        e_count[l][genome[i:i+1]] += 1

    state_2.reverse()

    i = 0
    while i < len(state_2):
        if state_2[i + 1] - state_2[i] < 50 - 1:
            del(state_2[i])
            del(state_2[i])
        else:
            i += 2

    e1_total = e1_count['A'] + e1_count['C'] + e1_count['G'] + e1_count['T']
    e2_total = e2_count['A'] + e2_count['C'] + e2_count['G'] + e2_count['T']

    for letter, count in e1_count.iteritems():
        e[1][letter] = 1.0 * e1_count[letter] / e1_total

    for letter, count in e2_count.iteritems():
        e[2][letter] = 1.0 * e2_count[letter] / e2_total

    a[1][1] = 1.0 * a1_count[1] / (a1_count[1] + a1_count[2])
    a[1][2] = 1.0 - a[1][1]
    a[2][1] = 1.0 * a2_count[1] / (a2_count[1] + a2_count[2])
    a[2][2] = 1.0 - a[2][1]

    return 1.0 * e1_total / (e1_total + e2_total), v_l_i_plus_1, state_2
