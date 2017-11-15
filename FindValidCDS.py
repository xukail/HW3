import re

def find_valid_CDS():
    filename = 'data/GCF_000091665.1_ASM9166v1_genomic.gbff'
    f = open(filename, 'r')

    valid_chars = set('0123456789CDS. \t\n')
    line = f.readline().upper()
    res = []

    while line is not None:
        if line.startswith('ORIGIN'):
            break
        if line.startswith('     CDS'):
            if not any((c not in valid_chars) for c in line):
                start_and_end = map(int, re.findall(r'\d+', line))
                res.append(start_and_end)

        line = f.readline().upper()

    f.close()
    # print res
    return res

if __name__ == "__main__":
    find_valid_CDS()
