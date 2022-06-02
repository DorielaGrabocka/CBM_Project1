
def DC3_recursive(t_genome, alphabet_size):

    SA12 = [] #list holding position of suffixes, where position is not divisible by 3
    for i in range(len(t_genome)):
        if i % 3 != 0:
            SA12 += [i]
    SA12 = bucket_sort(t_genome, alphabet_size, SA12, 2)
    SA12 = bucket_sort(t_genome, alphabet_size, SA12, 1)
    SA12 = bucket_sort(t_genome, alphabet_size, SA12)

    new_alpha = collect_alphabet(t_genome, SA12)
    if len(new_alpha) < len(SA12):
        u = build_u(t_genome, new_alpha)
        sa_u = DC3_recursive(u, len(new_alpha) + 2)
        m = len(sa_u) // 2
        SA12 = []
        for i in sa_u:
            if i != m:
                SA12 += [u_index(i, m)]

    if len(t_genome) % 3 == 1:

        SA3 = [len(t_genome) - 1]
    else:
        SA3 = []
    for i in SA12:
        if i % 3 == 1:
            SA3 += [i - 1]
    SA3 = bucket_sort(t_genome, alphabet_size, SA3)
    return merge(t_genome, SA12, SA3)

#This function converts the characters of the original genome, to their corresponding number in the alphabet
#because DC3 algorithm needs numbers to work.
def DC3(genome):
    map_base_to_index = {"$":0 , "A":1, "C":2, "G":3, "T":4 }
    t_genome = []
    for c in genome:
        t_genome += [map_base_to_index[c]]
    return DC3_recursive(t_genome, 5)

#Make sure that we are alwaays inside boundaries.
def safe_index(t_genome, index):
    if index >= len(t_genome):
        return 0
    else:
        return t_genome[index]

#Here we count the frequencies of each letter(now corresponding to a number)in the genome
def count_letters(t_genome, alphabet_size):
    counts = [0] * alphabet_size
    for c in t_genome:
        counts[c] += 1
    return counts

def cumulative_sum(counts):
    res = [0] * len(counts)
    sum = 0
    for i, k in enumerate(counts):
        res[i] = sum
        sum += k
    return res

def bucket_sort(t_genome, alphabet_size, index, offset = 0):
    sort_symbs = [] #list of letters at position offset of the given index
    for i in index:
        sort_symbs += [safe_index(t_genome, i + offset)]
    counts = count_letters(sort_symbs, alphabet_size)
    buckets = cumulative_sum(counts)
    out = [0] * len(index) #here we create the buckets holding the indexes to be sorted
    for i in index:
        bucket = safe_index(t_genome, i + offset)
        out[buckets[bucket]] = i
        buckets[bucket] += 1
    return out


def collect_alphabet(t_genome, index):
    alphabet = {}
    for i in index:
        triplet = (safe_index(t_genome, i), safe_index(t_genome, i + 1), safe_index(t_genome, i + 2))
        if triplet not in alphabet:
            alphabet[triplet] = len(alphabet) + 2 # saves size of the new triplet alphabet
    return alphabet

def build_u(t_genome, alphabet):
    u = []
    for i in range(1, len(t_genome), 3):
        u += [alphabet[(safe_index(t_genome, i), safe_index(t_genome, i + 1), safe_index(t_genome, i + 2))]]
    u += [1]
    for i in range(2, len(t_genome), 3):
        u += [alphabet[(safe_index(t_genome, i), safe_index(t_genome, i + 1), safe_index(t_genome, i + 2))]]
    return u

def u_index(i, m):
    if i < m:
        return 1 + 3 * i
    else:
        return 2 + 3 * (i - m - 1)

def merge(t_genome, SA12, SA3):
    ISA = {}
    for i in range(len(SA12)):
        ISA[SA12[i]] = i
    SA = []
    i, j = 0, 0
    while i < len(SA12) and j < len(SA3):
        if is_smaller(t_genome, SA12[i], SA3[j], ISA):
            SA.append(SA12[i])
            i += 1
        else:
            SA.append(SA3[j])
            j += 1
    SA.extend(SA12[i:])
    SA.extend(SA3[j:])
    return SA

def is_smaller(t_genome, i, j, ISA):
    g_at_i = safe_index(t_genome, i)
    g_at_j = safe_index(t_genome, j)
    if g_at_i < g_at_j:
        return True
    if g_at_i > g_at_j:
        return False
    if i % 3 != 0 and j % 3 != 0:
        return ISA[i] < ISA[j]
    return is_smaller(t_genome, i + 1, j + 1, ISA)

