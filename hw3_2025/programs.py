from pprint import pprint

#borrowed from the internet
def read_fasta(file_path):
    sequence = []
    with open(file_path, "r") as file:
        for line in file:
            if not line.startswith(">"):  # Ignore header lines
                sequence.append(line.strip())  # Collect sequence lines
    
    return "".join(sequence)  # Combine into a single string

#borrowed from internet
def parse_fasta_reads(file_path):
    """Reads a FASTA file and returns a list of separate read sequences."""
    reads = []
    with open(file_path, "r") as file:
        seq = ""
        for line in file:
            if line.startswith(">"):  # Header line
                if seq:  # Store the previous read
                    reads.append(seq)
                seq = ""  # Reset for new read
            else:
                seq += line.strip()  # Append sequence lines
        if seq:  # Add last read
            reads.append(seq)
    
    return reads  # Returns a list of individual read sequences


def get_kmers(data, k = 10):
    kmers = dict()
    res = set()
    for i in range(len(data) - k + 1):
        kmer = data[i:i + k]
        res.add(kmer)
        if kmer not in kmers:
            kmers[kmer] = 1
        else:
            kmers[kmer] += 1
    
    #return the top 10 kmers
    res = sorted(res, key = lambda x: kmers[x], reverse = True) #borrowed
    return res


def classify_reads(reads_genome, kmers_genome):
    
    #print("KMERS GENOME", kmers_genome)
    counts = dict()
    for read in reads_genome:
        
        #classify the reads in the reads genome      
        top_kmers = get_kmers(read)
        top_kmers = set(top_kmers)
        best_match = None
        max_matches = 0

        #print("KMERS GENOME", kmers_genome)
        for genome in kmers_genome:
            kmers = kmers_genome[genome] #the top kmers for the genome
            intersection = top_kmers.intersection(kmers)
            
            if len(intersection) > max_matches:
                max_matches = len(intersection)
                best_match = genome
        
        counts[read] = best_match #each read is mapped to its best match and how much it overlaps, used to compute relative frequency
    return counts





def main():
    # Example usage
    file_path = "data/genome1.fa"
    genome1 = read_fasta(file_path)

    file_path = "data/genome2.fa"
    genome2 = read_fasta(file_path)

    file_path = "data/genome3.fa"
    genome3 = read_fasta(file_path)

    file_path = "data/genome4.fa"
    genome4 = read_fasta(file_path)

    file_path = "data/reads.fa"
    read_genome = parse_fasta_reads(file_path)




    #2a-b
    top_kmers1 = get_kmers(genome1)
    top_kmers2 = get_kmers(genome2)
    top_kmers3 = get_kmers(genome3)
    top_kmers4 = get_kmers(genome4)
    
    # kmers_genomes = {'genome1': top_kmers1, 
    #                  'genome2': top_kmers2, 
    #                  'genome3': top_kmers3, 
    #                  'genome4': top_kmers4}

    
    # read_counts = classify_reads(read_genome, kmers_genomes)

    # #now we need to state the relative frequency of each genome

    # freqs = {'genome1': 0, 'genome2': 0, 'genome3': 0, 'genome4': 0}    
    
    # for read in read_counts:
    #     assert(read in read_counts)
    #     genome = read_counts[read]

    #     print(genome)
    #     if genome != None:
    #         freqs[genome] += 1
    
    # total_reads = 1000

    # print("Genome 1 Frequency: ", freqs['genome1']/total_reads)
    # print("Genome 2 Frequency: ", freqs['genome2']/total_reads)
    # print("Genome 3 Frequency: ", freqs['genome3']/total_reads)
    # print("Genome 4 Frequency: ", freqs['genome4']/total_reads)

    #2C - D
       #now we have all the 10mers
    #we need to find the top ten in each that aren't present in any of the others

    set_1 = set(top_kmers1)
    set_2 = set(top_kmers2)
    set_3 = set(top_kmers3)
    set_4 = set(top_kmers4)

    disc1 = []
    disc2 = []
    disc3 = []
    disc4 = []

    #get the top ones in genome 1 that are discriminative
    for kmer in top_kmers1:
        if kmer not in set_2 and kmer not in set_3 and kmer not in set_4:
            disc1.append(kmer)
    
    # disc1_2c = disc1[:5]

    #get the top ones in genome 2 that are discriminative
    for kmer in top_kmers2:
        if kmer not in set_1 and kmer not in set_3 and kmer not in set_4:
            disc2.append(kmer)
    # disc2_2c = disc2[:5]

    #get the top ones in genome 3 that are discriminative
    for kmer in top_kmers3:
        if kmer not in set_1 and kmer not in set_2 and kmer not in set_4:
            disc3.append(kmer)
    # disc3_2c = disc3[:5]

    #get the top ones in genome 4 that are discriminative
    for kmer in top_kmers4:
        if kmer not in set_1 and kmer not in set_2 and kmer not in set_3:
            disc4.append(kmer)
    # disc4_2c = disc4[:5]

    #classify the reads this time with the discriminative ones.

    print("DISCRIMINATIVE KMERS")
    kmers_genomes = {'genome1': disc1, 
                     'genome2': disc2, 
                     'genome3': disc3, 
                     'genome4': disc4}

    
    read_counts = classify_reads(read_genome, kmers_genomes)
    
    freqs = {'genome1': 0, 'genome2': 0, 'genome3': 0, 'genome4': 0}    
    
    for read in read_counts:
        genome = read_counts[read]

        if genome != None:
            freqs[genome] += 1
    
    total_reads = 1000
    print("Genome 1 Frequency: ", freqs['genome1']/total_reads)
    print("Genome 2 Frequency: ", freqs['genome2']/total_reads)
    print("Genome 3 Frequency: ", freqs['genome3']/total_reads)
    print("Genome 4 Frequency: ", freqs['genome4']/total_reads)

main()
