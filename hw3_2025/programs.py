#borrowed from the internet
def read_fasta(file_path):
    sequence = []
    with open(file_path, "r") as file:
        for line in file:
            if not line.startswith(">"):  # Ignore header lines
                sequence.append(line.strip())  # Collect sequence lines
    
    return "".join(sequence)  # Combine into a single string

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


def get_kmers(data, k = 4):
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
    res = sorted(res, key = lambda x: kmers[x], reverse = True)
    res = res[:10]
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
        
        counts[read] = (best_match, max_matches) #each read is mapped to its best match and how much it overlaps
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


    
    
    top_kmers1 = get_kmers(genome1)
    top_kmers2 = get_kmers(genome2)
    top_kmers3 = get_kmers(genome3)
    top_kmers4 = get_kmers(genome4)
    
    

    
    kmers_genomes = {'genome1': top_kmers1, 
                     'genome2': top_kmers2, 
                     'genome3': top_kmers3, 
                     'genome4': top_kmers4}

    
    read_counts = classify_reads(read_genome, kmers_genomes)
    
    #now we need to state the relative frequency of each genome

    freqs = {'genome1': 0, 'genome2': 0, 'genome3': 0, 'genome4': 0}    
    for read in read_counts:
        genome, similars = read_counts[read]
        freqs[genome] += 1
    
    total_reads = 1000

    print(freqs)
    print("Genome 1 Frequency: ", freqs['genome1']/total_reads)
    print("Genome 2 Frequency: ", freqs['genome2']/total_reads)
    print("Genome 3 Frequency: ", freqs['genome3']/total_reads)
    print("Genome 4 Frequency: ", freqs['genome4']/total_reads)




main()
