from LookForORF import look_for_orfs
from FindValidCDS import find_valid_CDS
import collections
import math
import matplotlib.pyplot as plt


def main():
    orfs, genome = look_for_orfs()
    CDS_list = find_valid_CDS()

    print '\nOutputs:\n'
    # print_1a(orfs)
    # print_1bc(orfs)
    # print_1d(CDS_list)
    
    numPatterns_Q, threeLettersMap_Q, fourLettersMap_Q, orf_map = build_model(orfs, genome, False, CDS_list)
    numPatterns_P, threeLettersMap_P, fourLettersMap_P, orf_map = build_model(orfs, genome, True, CDS_list)
    # print_1e([fourLettersMap_P, fourLettersMap_Q])

    # calculate Markov model score and make table of summarie for each orf
    for key, value in orf_map.iteritems():
        orf_str = genome[key - 1: key + value['length']]
        calculate_MMS(orf_str, value, [threeLettersMap_P, fourLettersMap_P], 
        [threeLettersMap_Q, fourLettersMap_Q], numPatterns_P, numPatterns_Q)

    # print_1f(orf_map)

    # draw plot for Q2
    axes = plt.gca()
    axes.set_xlim([0.0,1.0])
    axes.set_ylim([0.0,1.0])

    # split map into 1000 coordiantes according to orf's length
    coordinates = []
    for key, value in sorted(orf_map.iteritems(), key=lambda (k,v): (v['length'],k)):
        print key, ' ', value

    print '' #end line

def build_model(orfs, genome, p_or_q, CDS_list):
    threeLettersMap = {}
    fourLettersMap = {}
    orf_map = {}
    numPatterns = 0
    for i in range(len(orfs)):
        for j in range(len(orfs[i])):
            orf_len = 0
            if j is 0:
                orf_len =  orfs[i][j] - i
            else:
                orf_len = orfs[i][j] - orfs[i][j - 1]
            
            judge = True
            if p_or_q:
                judge = orf_len > 1400

                orf_content = {}
                orf_content['length'] = orf_len
                orf_content['status'] = False
                for cds in CDS_list:
                    if cds[1] == orfs[i][j]:
                        orf_content['status'] = True
                orf_map[orfs[i][j] - orf_len + 1] = orf_content

            else:
                judge = orf_len < 50

            if judge:
                for l in range(orfs[i][j] - orf_len, orfs[i][j] - 6, 1):
                    if genome[l : l + 3] not in threeLettersMap:
                        threeLettersMap[genome[l:l+3]] = 0
                    if genome[l : l + 4] not in fourLettersMap:
                        fourLettersMap[genome[l:l+4]] = 0
                    fourLettersMap[genome[l:l+4]] += 1
                    threeLettersMap[genome[l:l+3]] += 1
                    numPatterns += 1
    return numPatterns, threeLettersMap, fourLettersMap, orf_map

# calculate the Markov model score for a given ORF, P maps & Q maps,
# and the sum of total number of pattern occurences in such maps
def calculate_MMS(orf_str, orf_content, maps_P, maps_Q, numPatterns_P, numPatterns_Q):
    orf_content['MMS'] = 0
    if orf_content['length'] == 3:
        return 0
    log_p_val = math.log(1.0 * maps_P[0][orf_str[:3]] / numPatterns_P, 2)
    log_q_val = math.log(1.0 * maps_Q[0][orf_str[:3]] / numPatterns_Q, 2)
    for i in range(3, orf_content['length'] - 3): # exclude stop
        pattern_4 = orf_str[i - 3: i + 1]
        pattern_3 = orf_str[i - 3: i]
        if pattern_3 in maps_P[0]:
            log_p_val += math.log(1.0 * maps_P[1][pattern_4] / maps_P[0][pattern_3], 2)
        else:
            log_p_val += math.log(0.5, 2)
        if pattern_3 in maps_Q[0]:
            log_q_val += math.log(1.0 * maps_Q[1][pattern_4] / maps_Q[0][pattern_3], 2)
        else:
            log_q_val += math.log(0.5, 2)
    mms = log_p_val - log_q_val
    orf_content['MMS'] = mms
    return mms  

def print_1a(orfs):
    print '1a'
    print '\tThe number of ORFs for:'
    print '\tstart positions have indices congruent to 1 mod 3: ', len(orfs[0])
    print '\tstart positions have indices congruent to 2 mod 3: ', len(orfs[1])
    print '\tstart positions have indices congruent to 0 mod 3: ', len(orfs[2]), '\n'

def print_1bc(orfs):
    count50 = 0
    count1400 = 0
    for i in range(len(orfs)):
        for j in range(len(orfs[i])):
            orf_len = 0
            if j is 0:
                orf_len =  orfs[i][j] - i
            else:
                orf_len = orfs[i][j] - orfs[i][j - 1]
            if orf_len < 50:
                count50 += 1
            if orf_len > 1400:
                count1400 += 1

    print '1b'
    print '\tThe total number of ORFs of length less than 50: ', count50, '\n'
    print '1c'
    print '\tThe total number of ORFs of length greater than 1400: ', count1400, '\n'

def print_1d(CDS_list):
    print '1d'
    print '\tThe total number of "simple forward strand CDSs" you found in GenBank: ', len(CDS_list)
    print ''

def print_1e(map_list):
    letters = {'A' : 0, 'C' : 1, 'G' : 2, 'T' : 3}
    pAAxy_table = [[0, 0, 0, 0],[0, 0, 0, 0],[0, 0, 0, 0],[0, 0, 0, 0]]
    print '1e'
    for m in range(2):
        if m == 0:
            print '\tP(A | Axy):'
        else:
            print '\tQ(A | Axy):'
        print '\t       A\tC\tG\tT'

        for key, value in map_list[m].iteritems():
            if key[0] == 'A' and key[3] == 'A':
                pAAxy_table[letters[key[1]]][letters[key[2]]] = value
        for letter, i in collections.OrderedDict(sorted(letters.items())).iteritems():
            print '\t   ', letter, ' {}\t{}\t{}\t{}\t'.format(
            pAAxy_table[i][0], pAAxy_table[i][1], pAAxy_table[i][2], pAAxy_table[i][3])
        print ''

def print_1f(orf_map):
    print '1f'
    read_50 = True
    read_1400 = True
    less_than_50 = []
    greater_than_1400 = []
    for key in sorted(orf_map.iterkeys()):
        value = orf_map[key]
        summary = ''
        if value['length'] < 50 or value['length'] > 1400:
            summary = 'start coordinate: ' + str(key) + ', length: ' + str(value['length']) + ', log-base-2 Markov model score: ' + str(value['MMS']) + ', simple forward strand CDSs: ' + str(value['status'])
        if value['length'] < 50 and read_50:
            less_than_50.append(summary)
        if value['length'] > 1400 and read_1400:
            greater_than_1400.append(summary)
        if len(less_than_50) == 5:
            read_50 = False
        if len(greater_than_1400) == 5:
            read_1400 = False
        if not read_50 and not read_1400:
            break;
    print '\tLength less than 50:'
    for summary in less_than_50:
        print '\t', summary
        
    print '\n\tLength greater than 1400:'
    for summary in greater_than_1400:
        print '\t', summary
    
    print ''

if __name__ == "__main__":
    main()