import mmh3
from bitarray import bitarray
from Bio import SeqIO


class AbaumanniiBloomfilter:
    # Implementation of the Bloomfilter Project for Acinetobacter baumannii
    # Code partly from https://github.com/Phelimb/BIGSI

    clonetypes = 8
    hits_per_filter = [0] * clonetypes
    array_size = 98501457
    hashes = 7
    k = 20
    names = ['IC1', 'IC2', 'IC3', 'IC4', 'IC5', 'IC6', 'IC7', 'IC8']
    lookup_ct = [True, True, True, True, True, True, True, True]
    number_of_kmeres = 1

    def __init__(self):

        # creating matrix for filters
        #
        self.matrix = bitarray(self.clonetypes * self.array_size)

        # initializing all with 0
        self.matrix.setall(False)

    def print_matrix(self):
        #  prints matrix

        counter = 0
        for i in range(self.clonetypes):
            print(self.matrix[counter: self.array_size + counter])
            counter = counter + self.array_size

    # Getter

    def set_lookup(self, new):
        self.lookup_ct = new

    def get_hits_per_filter(self):
        return self.hits_per_filter

    def get_kmeres_per_sequence(self):
        # returns number of k-meres per file
        return self.number_of_kmeres

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

        # get the col for each hash value
        for i in range(self.hashes):

            # taking a col with number of elements = number of clonetypes

            col = []
            counter = 0  # counter for jumping from row to row
            for y in range(self.clonetypes):

                # getting col by getting specific matrix positions from hash-values
                col.append(self.matrix[positions[i] + counter])
                counter += self.array_size


            # checking row for 0s and 1s
            for j in range(len(col)):

                # if that cell is True und no False before
                if col[j] & hits[j] :
                    pass
                else:
                    hits[j] = False

        # Updating hits per filter
        for i in range(len(hits)):

            # If hit is True
            if hits[i] & self.lookup_ct[i]:
                # Update hit counter
                self.hits_per_filter[i] += 1

    def train(self, kmer, clonetype):
        # trains specific filter for a k-mer

        # getting hash Values
        positions = self.hash(kmer)
        # changing 0s to 1 in filter
        for i in range(len(positions)):
            # getting position of cell
            self.matrix[self.array_size * clonetype + positions[i]] = True

    def train_sequence(self, filepath, clonetype):
        # trains whole sequence into filter

        # for each sequence (in multi-FASTA file)
        for sequence in SeqIO.parse(filepath, "fasta"):

            # for each k-mere
            for i in range(len(sequence.seq) - self.k + 1):
                # trains k-mere into filter
                self.train(str(sequence.seq[i: i + self.k]), clonetype)

    def lookup_sequence(self, path):
        # uses lookup function for whole sequence

        # Counter of k-meres
        self.number_of_kmeres = 0
        self.hits_per_filter = [0] * self.clonetypes

        # for each sequence (in multi-FASTA file)
        for sequence in SeqIO.parse(path, "fasta"):

            # Updating number of K-meres
            #self.number_of_kmeres += len(sequence.seq) - self.k + 1

            # for each k-mere
            for i in range(0, len(sequence.seq) - self.k + 1, 2):

                # lookup for all k-meres in filter
                self.lookup(str(sequence.seq[i: i + self.k]))
                self.number_of_kmeres += 1


    def add_filter(self, name):
        # adds one EMPTY(!) filter to the current filter
        # the new filter will be added on the back
        # use the train(seq) function for training

        # the new filter is a row that has to be inserted

        # don't forget to save the filter after adding/training it

        # creating new filter
        new = bitarray(self.array_size)

        # setting all bits to False
        new.setall(False)

        # adding filter to matrix
        self.matrix.extend(new)

        # updating ICs
        self.clonetypes += 1

        # appending new name of filter
        self.names.append(name)

        # adding True to lookup
        self.lookup_ct.append(True)


    def remove_filter(self, ct):
        # deletes filter
        # !!! Counting starts with 0
        # also deletes name from name-list

        del self.names[ct]
        del self.lookup_ct[ct]
        self.clonetypes -= 1

        # always deleting same position with pop()
        #

        for i in range(self.array_size):
            self.matrix.pop(ct * self.array_size)

    def get_score(self):
        # calculates score for all clonetypes
        # Score is #hits / #kmeres

        score = []

        # calculates float for each value in [hits per filter]
        for i in range(self.clonetypes):
            score.append(round(float(self.hits_per_filter[i]) / float(self.number_of_kmeres), 2))


        return score

    def get_names(self):
        return self.names

    def lookup_fastq(self, path, number_of_reads):
        # lookup for fastq files,
        # Source: http://biopython.org/DIST/docs/tutorial/Tutorial.html , 20.1.11  Indexing a FASTQ file
        # read fastq file

        fq_dict = SeqIO.index(path, "fastq")
        # getting number of reads in file
        reads = len(fq_dict)

        # counter for kmeres
        self.number_of_kmeres = 0
        seqs = len(fq_dict)
        c = 0
        # lookup for all reads
        for i in fq_dict.keys():
            if c == number_of_reads:
                break
            c += 1
            # getting read
            single_read = fq_dict[i].seq
            # lookup for k-meres in read
            for j in range(len(single_read)):
                # updating counter
                self.number_of_kmeres += 1

                # lookup for kmer
                self.lookup(str(single_read[j: j + self.k]))

