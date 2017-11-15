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

    res = look_for_stops(genome)
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

if __name__ == "__main__":
    main()