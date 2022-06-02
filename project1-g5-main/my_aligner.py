import pickle
import numpy as np
import math
import sys

from aligner import *
from DC3 import *

class FMIndexer(ReferenceIndex):
    # __________________our get_seed
    ''' A class for indexing the reference genome to allow for fast searches. '''

    bwt = None
    C = None
    occ = None
    index = None
    seed_length = None

    def __init__(self, ref_file: str) -> None:
        self.ref_file = ref_file + ".fa"
        self.fq1 = ref_file + ".sim1.fq"
        self.fq2 = ref_file + ".sim2.fq"
        self.ref_concat = ref_file + ".sim_errFree.fa"
        self.genome = self.get_genome_string_from_file(ref_file)
        self.length_of_original_genome = len(self.genome)
        self.n_matrix, self.genome = self.compute_n_matrix_and_substitute_n(self.genome)
        super().__init__(ref_file)

    def get_genome_string_from_file(self, ref_file) -> str:
        genome = ''
        with open(ref_file, 'r') as rf:
            for l in rf:
                if l[0] != '>':
                    genome += l.rstrip()
        return genome

    # This method is used to generate the index of the reference genome
    def _generate_index(self) -> None:

        self.index = DC3(self.genome+'$')
        self.bwt = self.find_Bwt_letters(self.index, self.genome)
        self.occ = self.build_occ_matrix(
            "".join(self.bwt))  # joined with "" to convert the array of characters to string
        self.C = self.build_c_matrix(self.occ)
        self.serialize(self.ref_file.filename.decode() + ".idx")

    # we iterate once and build n matrix and remove Ns
    # Everything in 1xO(n) time
    def compute_n_matrix_and_substitute_n(self, genome):
        array_of_ns = []  # array of key-value pairs, the key s the position, the value is the number of N's starting at that position
        length = len(self.genome)
        i = 0
        original_index = 0
        while i < length:
            foundN = False
            counter = 0
            while i < length and self.genome[i] == 'N':
                counter = counter + 1
                foundN = True
                i = i + 1

            if (counter != 0):
                array_of_ns.append((original_index, counter))
                self.genome = self.genome[:i - counter] + "$" + self.genome[i:]
                length = length - counter + 1

            if foundN:
                original_index = original_index + counter
            else:
                original_index = original_index + 1
            i = i - counter + 1

        return array_of_ns, self.genome

    def calculate_real_position_in_original(self, n_matrix, index, length_of_genome):
        sum = 0  # holds number of N up to a certain position
        pos = 0  # holds the pointer position taking into account $
        no_of_iterations = len(n_matrix)
        if len(n_matrix) == 0: return index
        for i in range(0, no_of_iterations):
            a = n_matrix[i][0]
            b = n_matrix[i][1]
            if (i + 1) == no_of_iterations:
                next_index = length_of_genome - 1
            else:
                next_index = n_matrix[i + 1][0]
            n = next_index - (a + b)  # here we count the number of non-N characters between two groups of N
            pos = pos + n + 1
            sum = sum + b
            if pos > index:
                return index + sum - (i + 1)  # here is the original index

    # Function used to compute and return all suffixes of a genome
    def get_suffixes(self, word):
        suffixes = []
        word += '$'
        for i in range(len(word)):
            suffixes.append(word[len(word) - i - 1:])

        return suffixes

    # Function used to do the suffix array
    def compute_suffix_array(self, word):
        all_suffixes = self.get_suffixes(word)
        sorted_suffixes = sorted(all_suffixes)
        suffix_index_array = []
        for el in sorted_suffixes:
            suffix_index_array.append(len(all_suffixes) - all_suffixes.index(el) - 1)
        return suffix_index_array

    # function used to find the BWT transformation of the genome string
    def find_Bwt_letters(self, suffix_array, word):
        letters = []
        word = word + '$'
        for index in suffix_array:
            letters.append(word[index % len(word) - 1])
        return letters

    # Function used to build the occ matrix
    def build_occ_matrix(self, genome_string):
        # $->0, A->1, C->2, G->3, T->4
        column_count = len(genome_string)
        row_count = 5
        character_array = {'$': 0, 'A': 0, 'C': 0, 'G': 0, 'T': 0}
        occ_matrix = [[0] * column_count for _ in range(row_count)]
        for i in range(column_count):
            character_array[genome_string[i]] = character_array[genome_string[i]] + 1
            occ_matrix[0][i] = character_array['$']
            occ_matrix[1][i] = character_array['A']
            occ_matrix[2][i] = character_array['C']
            occ_matrix[3][i] = character_array['G']
            occ_matrix[4][i] = character_array['T']
        return occ_matrix

    # Function used to build the count matrix
    def build_c_matrix(self, occ_matrix):
        cols = len(occ_matrix[0]) - 1
        # order: $,A,C,G,T
        c_matrix = [0,
                    occ_matrix[0][cols],
                    occ_matrix[1][cols] + occ_matrix[0][cols],
                    occ_matrix[2][cols] + occ_matrix[1][cols] + occ_matrix[0][cols],
                    occ_matrix[3][cols] + occ_matrix[2][cols] + occ_matrix[1][cols] + occ_matrix[0][cols]]
        return c_matrix

    # The seed is the pattern we are going to match.
    # The algorithm will return the left and the right boundaries were a match might be located.
    def find_hit_range(self, query_seed):
        map_base_to_index = {"$": 0, "A": 1, "C": 2, "G": 3, "T": 4}
        i = len(query_seed) - 1  # i set to the length of the seed
        c = query_seed[i]  # initialize c to the last charcter of a seed
        sp = self.C[map_base_to_index[c]]  # left boundary where the search will start(of ref_index)
        # right boundary(of ref index)
        if c != "T":
            ep = self.C[(map_base_to_index[c] + 1) % 5] - 1
        else:
            ep = len(self.genome)
        while sp <= ep and i >= 1:
            c = query_seed[i - 1]
            sp = self.C[map_base_to_index[c]] + self.occ[map_base_to_index[c]][sp - 1]
            ep = self.C[map_base_to_index[c]] + self.occ[map_base_to_index[c]][ep] - 1
            i = i - 1
        return sp, ep

    # This method is used to  write the index  in a file
    def serialize(self, out_file: str) -> None:
        sys.stderr.write("Writing index to {}\n".format(out_file))
        with open(out_file, "wb") as f:
            pickle.dump(self.index, f)
            pickle.dump(self.bwt, f)
            pickle.dump(self.occ, f)
            pickle.dump(self.C, f)

    # This method is used to load the index from the file
    def load(self, in_file: str) -> None:
        #assert(False)
        with open(in_file, "rb") as f:
            self.ref_file = pysam.FastaFile(".".join(in_file.split(".")[:-1]))
            sys.stderr.write("Loading index from {}\n".format(in_file))
            self.index = pickle.load(f)
            self.bwt = pickle.load(f)
            self.occ = pickle.load(f)
            self.C = pickle.load(f)


    # !!!!!!!!!!!!!!!!!NEW get_seeds_method here!!!!!!!!!!!!!!!!!!!!!!!####
    def get_seeds(self, query_record: SeqRecord):
        seeds = []
        seeds_reversed = []  # list of seeds that we found for this read
        complete_read = query_record.seq  # original read
        complete_read_reversed = query_record.reverse_complement().seq  # reversed read
        length_of_read = len(query_record.seq)
        number_of_levels = math.floor(math.log(length_of_read, 2))
        # number of levels must be included
        for i in range(0, number_of_levels + 1):
            power_of_two = int(math.pow(2, i))
            offset = int(length_of_read / power_of_two) - 1
            start_pos = 0
            end_pos = start_pos + offset
            no_of_children = power_of_two
            for j in range(no_of_children):
                seed_query = complete_read[start_pos:end_pos + 1]  # because last index is not inclusive
                seed_query_reversed = complete_read_reversed[start_pos:end_pos + 1]
                ref_start, ref_end = self.find_hit_range(seed_query)
                ref_start_rev, ref_end_rev = self.find_hit_range(seed_query_reversed)
                ref_end = min(ref_end, ref_start + 3)
                ref_end_rev = min(ref_end_rev, ref_start_rev + 3)
                for k in range(ref_start, ref_end + 1):
                    index_in_original = self.calculate_real_position_in_original(self.n_matrix,self.index[k], self.length_of_original_genome)
                    seed_object = Seed(query_record.id,
                                       seed_query, False,
                                       start_pos, end_pos + 1,
                                       index_in_original,
                                       index_in_original + len(seed_query))
                    seeds.append(seed_object)

                for k in range(ref_start_rev, ref_end_rev + 1):
                    index_in_original = self.calculate_real_position_in_original(self.n_matrix, self.index[k], self.length_of_original_genome)
                    seed_object_reversed = Seed(query_record.id,
                                                seed_query_reversed, True,
                                                start_pos, end_pos + 1,
                                                index_in_original,
                                                index_in_original + len(seed_query_reversed))
                    seeds_reversed.append(seed_object_reversed)

                start_pos = end_pos + 1
                end_pos = start_pos + offset if j + 2 < no_of_children else length_of_read - 1
            # here we found the longest matching seed
            if len(seeds_reversed) > 0 or len(seeds) > 0:
                break

        if len(seeds) > len(seeds_reversed) and len(seeds_reversed) != 0:
            return seeds
        elif len(seeds) < len(seeds_reversed) and len(seeds) != 0:
            return seeds_reversed
        else:
            seeds_reversed.extend(seeds)
            return seeds_reversed


class MyExtender(Extender):
    buffer = 2
    match = 2
    mismatch = -1
    gap = -3
    bound = 7

    def extend(self, index: ReferenceIndex, query_record: SeqRecord) -> Alignment:
        assert (not index.ref_file is None)
        best_alignment = None
        best_align_score = -1000000
        for seed in index.get_seeds(query_record):
            assert (not seed.orientation is None)
            assert (not seed.ref_start is None)
            assert (not seed.ref_end is None)
            assert (not seed.query_start is None)
            assert (not seed.query_end is None)
            query = query_record.seq if not seed.orientation else query_record.reverse_complement().seq
            ref_start = max(seed.ref_start - seed.query_start - self.buffer, 0)
            ref_end = min(seed.ref_end + len(query) - seed.query_end + self.buffer, index.ref_file.lengths[0])
            ref = index.fetch_ref(ref_start, ref_end)
            alignment = self.align(seed, ref, ref_start, ref_end, query)
            if alignment.score > best_align_score:
                best_align_score = alignment.score
                best_alignment = alignment

        return best_alignment

    def align(self, seed, ref, ref_start, ref_end, query):
        score = self.get_score(query, ref)
        alignment = Alignment(seed.query_id, seed.query_seq, seed.orientation, seed.query_start, seed.query_end,
                              ref_start + self.buffer, ref_end - self.buffer, score, '')
        return alignment

    def get_score(self, query_record, reference_seq):
        ref_size = len(reference_seq) + 1
        query_size = len(query_record) + 1
        S = np.full((ref_size, query_size), -10000, dtype=int)
        for j in range(self.bound + 1):
            S[0, j] = self.gap * j
        for i in range(1, ref_size):
            if i <= self.bound:
                S[i, 0] = 0
            lb = max(1, i - self.bound)
            ub = min(query_size, i + 1 + self.bound)
            for j in range(lb, ub):
                S[i, j] = max(S[i - 1, j] + self.gap, S[i, j - 1] + self.gap,
                              S[i - 1, j - 1] + self.diff(reference_seq[i - 1], query_record[j - 1]))

        return max(S[:, query_size - 1])

    def diff(self, a, b):
        if a == b:
            return self.match
        else:
            return self.mismatch
