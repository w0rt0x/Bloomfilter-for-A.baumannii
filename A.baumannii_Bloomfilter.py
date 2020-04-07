import mmh3
from bitarray import bitarray
from Bio import SeqIO

class AbaumanniiBloomfilter:
    # Implementation of the Bloomfilter Project for Acinetobacter baumannii
    # Code partly from https://github.com/Phelimb/BIGSI

    clonetypes = 8
    hits_per_filter = [0] * clonetypes
    array_size = 0
    set_size = 0
    hashes = 0
    false_positive_rate = 0
    k = 20

    def __init__(self):

        # creating matrix for filters
        #
        self.matrix = bitarray(self.clonetypes * self.array_size)

        # initializing all with 0
        self.matrix.setall(False)

    # Getter

    def get_hits_per_filter(self):
        return self.hits_per_filter

    def get_kmeres_per_sequence(self):
        # returns number of k-meres per file
        return self.number_of_kmeres

    # Setter

    def set_k(self, new):
        # sets k-mer size
        self.k = new

    def reset_hits(self):
        self.hits_per_filter = [0] * self.clonetypes

    def set_matrix(self):
        # re-creates matrix

        self.matrix = bitarray(self.array_size * self.clonetypes)
        self.matrix.setall(False)

    def set_clonetypes(self, new):
        # Setter for number of clonetypes

        # and Hits per filter
        self.clonetypes = new
        self.hits_per_filter = [0] * self.clonetypes

        # re-creating matrix because parameter changed
        self.set_matrix()

    def set_array_size(self, new):
        # Setter for Array Size

        self.array_size = new

        # re-creating matrix because parameter changed
        self.set_matrix()

    def set_set_size(self, new):
        # Setter for Set Size

        self.set_size = new

    def set_hashes(self, new):
        # Setter for number of hash functions

        self.hashes = new

    # File management

    def save_clonetypes(self, path):
        # saving filters of clonetypes

        # creating file and saving matrix with the bitarray modul
        with open(path, 'wb') as fh:

            # writing to file with bitarray command
            self.matrix.tofile(fh)

    def read_clonetypes(self, path):
        # reading filters from file

        # Source:
        # https://stackoverflow.com/questions/6266330/python-bitarray-to-and-from-file

        temp = bitarray()

        with open(path, 'rb') as fh:
            temp.fromfile(fh)
        self.matrix = temp

    # Bloomfilter

    def hash(self, kmer):
        # Hashes given string and returns Positions for the Array

        # Empty list for Array positions
        positions = []

        # Creating hashes for needed number of hash functions
        for i in range(self.hashes):
            # mmh3 takes that string and a seed,
            # each hash function takes an individual seed
            # after that, the hash-value will me divided by the array size until
            # a position in the array is guaranteed
            positions.append(mmh3.hash(kmer, i) % self.array_size)

        return positions

    def lookup(self, kmer):
        # checks if an element is in the filters, returns list with True/False
        # takes kmer input string and checks all clonetypes if the k-mer is inside that set of kmers

        # getting positions
        positions = self.hash(kmer)

        # control if element is in filter
        hits = [True] * self.clonetypes


        # get the row for each hash value
        for i in range(self.hashes):

            # taking a row with number of elements = number of clonetypes
            row = self.matrix[positions[i] * self.clonetypes: positions[i] * self.clonetypes + self.clonetypes]

            # checking row for 0s and 1s
            for j in range(len(row)):

                # if that cell is True und no False before
                if row[j] & hits[j]:
                    pass
                else:
                    hits[j] = False

        # Updating hits per filter
        for i in range(len(hits)):

            # If hit is True
            if hits[i]:

                # Update hit counter
                self.hits_per_filter[i] += 1

    def train(self, kmer, clonetype):
        # trains specific filter for a k-mer

        # getting hash Values
        positions = self.hash(kmer)

        # changing 0s to 1 in filter
        for i in range(len(positions)):

            # getting position of cell
            self.matrix[positions[i] * self.clonetypes + clonetype] = True

    def train_sequence(self, filepath, clonetype):
        # trains whole sequence into filter

        # for each sequence (in multi-FASTA file)
        for sequence in SeqIO.parse(filepath, "fasta"):

            # for each k-mere
            for i in range(len(sequence.seq) - self.k + 1):

                # trains k-mere into filter
                self.train(str(sequence.seq[i : i + self.k]), clonetype)
                

    def lookup_sequence(self, path):
        # uses lookup function for whole sequence

        # Counter of k-meres
        self.number_of_kmeres = 0
        # for each sequence (in multi-FASTA file)
        for sequence in SeqIO.parse(path, "fasta"):

            # Updating number of K-meres
            self.number_of_kmeres += len(sequence.seq) - self.k + 1

            # for each k-mere
            for i in range(len(sequence.seq) - self.k + 1):

                # trains k-mere into filter
                self.lookup(str(sequence.seq[i : i + self.k]))
   
    def lookup_fastq(self, path):
        # lookup for fastq files,
        # Source: http://biopython.org/DIST/docs/tutorial/Tutorial.html , 20.1.11  Indexing a FASTQ file

        # read fastq file
        fq_dict = SeqIO.index(path, "fastq")
        # getting number of reads in file
        reads = len(fq_dict)

        # counter for kmeres
        self.number_of_kmeres = 0
        seqs = len(fq_dict)
        print('Reads: ' , seqs)
        c = 0                                                                                                                   
	# lookup for just 100.000 reads                                                                                                  
	for i in fq_dict.keys():
            if c == 100000:
                break
            c += 1
            # getting read
            single_read = fq_dict[i].seq
            #print(single_read)
            #lokkup for k-meres in read
            for j in range(len(single_read)):

                # updating counter
                self.number_of_kmeres += 1

                # lookup for kmer
                self.lookup(str(single_read[j: j + self.k]))
