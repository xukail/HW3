from LookForORF import look_for_orfs
from FindValidCDS import find_valid_CDS
# from sklearn.metrics import roc_curve, auc
import collections
import math
import matplotlib.pyplot as plt

def question5_2():
    orfs, genome = look_for_orfs()
    CDS_list = find_valid_CDS()
    plt.figure()
    
    numPatterns_Q, threeLettersMap_Q, fourLettersMap_Q, orf_map = build_model(orfs, genome, False, CDS_list)
    numPatterns_P, threeLettersMap_P, fourLettersMap_P, orf_map = build_model(orfs, genome, True, CDS_list)
    # print_1e([fourLettersMap_P, fourLettersMap_Q])

    # calculate Markov model score and make table of summarie for each orf

    num_true = 0
    A_count = 0
    B_count = 0
    x_true = []
    y_true = []
    x_false = []
    y_false = []
    x_A = []
    x_B = []
    y_A = []
    y_B = []
    for key, value in orf_map.iteritems():
        orf_str = genome[key - 1: key + value['length']]
        calculate_MMS(orf_str, value, [threeLettersMap_P, fourLettersMap_P], 
        [threeLettersMap_Q, fourLettersMap_Q], numPatterns_P, numPatterns_Q)
        # scatter points
        if value['status']:
            num_true += 1
            if value['length'] > 1400 or value['length'] < 50:
                x_true.append(value['length'])
                y_true.append(value['MMS'])
        else:
            if value['length'] > 1400 or value['length'] < 50:
                x_false.append(value['length'])
                y_false.append(value['MMS'])
        if value['length'] > 1400:
            x_B.append(value['length'])
            y_B.append(value['MMS'])
        elif value['length'] < 50:
            x_A.append(value['length'])
            y_A.append(value['MMS'])
    num_false = len(orf_map) - num_true

    x_A, y_A, x_y_A = sort_coordiante(x_A, y_A)
    x_B, y_B, x_y_B = sort_coordiante(x_B, y_B)
    
    A_median = [x_A[len(x_A) / 2], sorted(y_A)[len(y_A) / 2]]
    B_median = [x_B[len(x_B) / 2], sorted(y_B)[len(y_B) / 2]]

    # calculate the slope and y-intercept of cut line
    cut_x = 0.2 * (B_median[0] - A_median[0]) + A_median[0]
    slope0 = 1.0 * (B_median[1] - A_median[1]) / (B_median[0] - A_median[0])
    slope1 = -1.0 / slope0
    b0 = A_median[1] - slope0 * A_median[0]
    b1 = cut_x * (slope0 - slope1) + b0

    
    min_threshold_80_tpr_cut = {'b': 24350}
    fpr_list = []
    tpr_list = []
    b = -150
    calculate = True
    while calculate:
        count_fp = 0
        count_tp = 0
        for key, value in orf_map.iteritems():
            if value['MMS'] >= slope1 * value['length'] + b: # Predict True
                if value['status']:
                    count_tp += 1
                else:
                    count_fp += 1
        tpr = 1.0 * count_tp / num_true
        fpr = 1.0 * count_fp / num_false
        fpr_list.append(fpr)
        tpr_list.append(tpr)
        if tpr >= 0.8:
            min_threshold_80_tpr_cut['b'] = b
            min_threshold_80_tpr_cut['tp_count'] = count_tp
            min_threshold_80_tpr_cut['fp_count'] = count_fp
            min_threshold_80_tpr_cut['tpr'] = tpr
        if tpr == 0.0 and fpr == 0.0:
            calculate = False
        if b == 10650:
            print count_fp
        if tpr < 0.004:
            b += 10000
        elif tpr < 0.05:
            b += 1500
        elif tpr < 0.1:
            b += 500
        else:
            b += 100
    # calculate AUC for the ROC
    area_cut = calculate_area(fpr_list, tpr_list)

    print min_threshold_80_tpr_cut
    plt.figure()
    lw = 2
    plt.plot(fpr_list, tpr_list, color='red', lw=lw, label='Cut line shifting (AUC = %0.2f)' % area_cut)
    plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('ROC curves for Question 5')
    plt.legend(loc="lower right")
    plt.show()

def question5_1():
    orfs, genome = look_for_orfs()
    CDS_list = find_valid_CDS()
    plt.figure()
    
    numPatterns_Q, threeLettersMap_Q, fourLettersMap_Q, orf_map = build_model(orfs, genome, False, CDS_list)
    numPatterns_P, threeLettersMap_P, fourLettersMap_P, orf_map = build_model(orfs, genome, True, CDS_list)
    # print_1e([fourLettersMap_P, fourLettersMap_Q])

    # calculate Markov model score and make table of summarie for each orf

    num_true = 0
    A_count = 0
    B_count = 0
    x_true = []
    y_true = []
    x_false = []
    y_false = []
    x_A = []
    x_B = []
    y_A = []
    y_B = []
    for key, value in orf_map.iteritems():
        orf_str = genome[key - 1: key + value['length']]
        calculate_MMS(orf_str, value, [threeLettersMap_P, fourLettersMap_P], 
        [threeLettersMap_Q, fourLettersMap_Q], numPatterns_P, numPatterns_Q)
        # scatter points
        if value['status']:
            num_true += 1
            # if value['length'] > 1400 or value['length'] < 50:
            x_true.append(value['length'])
            y_true.append(value['MMS'])
        else:
            # if value['length'] > 1400 or value['length'] < 50:
            x_false.append(value['length'])
            y_false.append(value['MMS'])
        if value['length'] > 1400:
            x_B.append(value['length'])
            y_B.append(value['MMS'])
        elif value['length'] < 50:
            x_A.append(value['length'])
            y_A.append(value['MMS'])
    num_false = len(orf_map) - num_true

    x_A, y_A, x_y_A = sort_coordiante(x_A, y_A)
    x_B, y_B, x_y_B = sort_coordiante(x_B, y_B)
    
    A_median = [x_A[len(x_A) / 2], sorted(y_A)[len(y_A) / 2]]
    B_median = [x_B[len(x_B) / 2], sorted(y_B)[len(y_B) / 2]]

    # calculate the slope and y-intercept of cut line
    cut_x = 0.2 * (B_median[0] - A_median[0]) + A_median[0]
    slope0 = 1.0 * (B_median[1] - A_median[1]) / (B_median[0] - A_median[0])
    slope1 = -1.0 / slope0
    b0 = A_median[1] - slope0 * A_median[0]
    b1 = cut_x * (slope0 - slope1) + b0

    AB_x = [A_median[0], B_median[0]]
    AB_y = [A_median[1], B_median[1]]

    print slope1
    print b1
    print num_false

    plt.scatter([x_true], y_true, color='orange', s=5, label='True protein ORFs')
    plt.scatter(x_false, y_false, color='blue', s=5, label='Non-protein ORFs')

    plt.scatter([A_median[0]], [A_median[1]], color='red', s=30, label='A')
    plt.scatter([B_median[0]], [B_median[1]], color='green', s=30, label='B')

    plt.plot(AB_x, AB_y, color='black', lw=2, label='AB line segment')

    temp_x = [-1000.0, 10000.0]
    temp_y = [(-1000.0 * slope1 + b1) , (10000.0 * slope1 + b1)]
    plt.plot(temp_x, temp_y, color='yellow', lw=2, label='Cut AB line segment')
    plt.plot([50, 50], [-3000, 3000], color='grey', lw=2, label='x = 50')
    plt.plot([1400, 1400], [-3000, 3000], color='grey', lw=2, label='x = 1400')

    # fpr_list, tpr_list, min_threshold_80_tpr = get_fpr_and_tpr_list(orf_map, length_occur, mms_occur, num_true, num_false)


    lw = 2
    # plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
    plt.xlim([0.0, 3000])
    plt.ylim([0.00, 3000])
    plt.xlabel('Length')
    plt.ylabel('Markov Model Score')
    plt.title('Question 5')
    plt.legend(loc="upper right")
    plt.show()

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
    length_occur = []
    mms_occur = []
    num_true = 0
    for key, value in orf_map.iteritems():
        orf_str = genome[key - 1: key + value['length']]
        calculate_MMS(orf_str, value, [threeLettersMap_P, fourLettersMap_P], 
        [threeLettersMap_Q, fourLettersMap_Q], numPatterns_P, numPatterns_Q)
        if value['status'] == True:
            num_true += 1
        if value['length'] not in length_occur:
            length_occur.append(value['length'])
    num_false = len(orf_map) - num_true
    length_occur = sorted(length_occur)

    for key, value in sorted(orf_map.iteritems(), key=lambda (k,v): (v['MMS'],k)):
        if len(mms_occur) == 0:
            mms_occur.append(value['MMS'])
        elif value['MMS'] > mms_occur[len(mms_occur) - 1]:
            mms_occur.append(value['MMS'])

    # print_1f(orf_map)

    # draw plot for Q2
    # split map into 1000 coordiantes according to orf's length
    fpr_list, tpr_list, min_threshold_80_tpr = get_fpr_and_tpr_list(orf_map, length_occur, mms_occur, num_true, num_false)

    area_length = calculate_area(fpr_list[0], tpr_list[0])
    area_MMS = calculate_area(fpr_list[1], tpr_list[1])

    plt.figure()
    lw = 2
    plt.plot(fpr_list[0], tpr_list[0], color='red', lw=lw, label='Length (AUC = %0.2f)' % area_length)
    plt.plot(fpr_list[1], tpr_list[1], color='green', lw=lw, label='Markov Model Score (AUC = %0.2f)' % area_MMS)
    plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('ROC curves for Question 2')
    plt.legend(loc="lower right")
    plt.show()

    # print result for Q3 and Q4
    print min_threshold_80_tpr[0]
    print min_threshold_80_tpr[1]


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

def get_fpr_and_tpr_list(orf_map, length_occur, mms_occur, num_true, num_false):
    min_threshold_80_tpr_length = {'length': 8000}
    min_threshold_80_tpr_MMS = {'MMS': 100.0}
    fpr_list = [[],[]]
    tpr_list = [[],[]]
    for length_filter in length_occur:
        count_fp = 0
        count_tp = 0
        for key, value in orf_map.iteritems():
            if value['length'] >= length_filter: # Predict True
                if value['status']:
                    count_tp += 1
                else:
                    count_fp += 1
        tpr = 1.0 * count_tp / num_true
        fpr = 1.0 * count_fp / num_false
        fpr_list[0].append(fpr)
        tpr_list[0].append(tpr)
        if value['length'] <= min_threshold_80_tpr_length['length'] and tpr >= 0.8:
            min_threshold_80_tpr_length['length'] = value['length']
            min_threshold_80_tpr_length['tp_count'] = count_tp
            min_threshold_80_tpr_length['fp_count'] = count_fp
            min_threshold_80_tpr_length['tpr'] = tpr

    step = 80
    i = 0
    while i < len(mms_occur):
        count_fp = 0
        count_tp = 0
        for key, value in orf_map.iteritems():
            if value['MMS'] >= mms_occur[i]: # Predict True
                if value['status']:
                    count_tp += 1
                else:
                    count_fp += 1
        tpr = 1.0 * count_tp / num_true
        fpr = 1.0 * count_fp / num_false
        fpr_list[1].append(fpr)
        tpr_list[1].append(tpr)
        if value['MMS'] <= min_threshold_80_tpr_MMS['MMS'] and tpr >= 0.8:
            min_threshold_80_tpr_MMS['MMS'] = value['MMS']
            min_threshold_80_tpr_MMS['tp_count'] = count_tp
            min_threshold_80_tpr_MMS['fp_count'] = count_fp
            min_threshold_80_tpr_MMS['tpr'] = tpr
        if mms_occur[i] > -6:
            step = 10
        i += step

    return fpr_list, tpr_list, [min_threshold_80_tpr_length, min_threshold_80_tpr_MMS]

def calculate_area(x_list, y_list):
    x_to_y = {}
    for i in range(len(x_list)):
        x_to_y[x_list[i]] = y_list[i]
    area = 0.0
    x = []
    y = []
    for key in sorted(x_to_y.iterkeys()):
        x.append(key)
        y.append(x_to_y[key])
    for i in range(1, len(x)):
        area += (y[i] + y[i - 1]) * (x[i] - x[i - 1]) / 2.0
    return area

def sort_coordiante(x, y):
    x_y = []
    for i in range(len(x)):
        x_y.append([x[i], y[i]])

    x_y = sorted(x_y, key = lambda c: int(c[0]))
    res_x = []
    res_y = []
    for xy in x_y:
        res_x.append(xy[0])
        res_y.append(xy[1])
    return res_x, res_y, x_y

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
    # main()
    question5_1()
    # question5_2()
    print 'Done.'