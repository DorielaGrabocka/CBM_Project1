from re import A
import numpy as np
from aligner import *
import tqdm  # type: ignore
from Bio import pairwise2  # type: ignore
from Bio.Seq import Seq
import pysam

def score(ref: str, query: str) -> float:
    match = 0
    mismatch = -1
    gap_open = -1
    gap_extend = -1
    aln = pairwise2.align.globalmd(ref, query, match, mismatch,
                                   gap_open, gap_extend,
                                   gap_open, gap_extend,
                                   penalize_end_gaps=(True, False),
                                   one_alignment_only=True)
    fwd = -aln[0].score if len(aln) else float(len(ref))
    rc = Seq(query)
    aln = pairwise2.align.globalmd(ref, str(rc.reverse_complement()), match, mismatch,
                                   gap_open, gap_extend,
                                   gap_open, gap_extend,
                                   penalize_end_gaps=(True, False),
                                   one_alignment_only=True)
    rev = -aln[0].score if len(aln) else float(len(ref))
    
    return min(fwd, rev)

def eval(make_index: ReferenceIndex, make_extender: Extender, *args) -> None:
    """
    To measure performance with a ground truth
    :param make_index:
    :param make_extender:
    :param args:
    :return:
    """
 
    ref_prefixes = [
        # "data/chr22",
        #"data_small/genome.chr22.5K",
        "data_big/genome.chr22"
    ]
 
    for ref_prefix in ref_prefixes:
        ref = ref_prefix + ".fa"
        f1 = "data_big/output_5xCov1.fq"
        f2 = "data_big/output_5xCov2.fq"
        #f1 = ref_prefix +".sim1.fq"
        #f2 = ref_prefix +".sim2.fq"
        #sam_file = ref_prefix +".sim_errFree.sam"
        #sam = pysam.AlignmentFile(sam_file, 'r')
 
        index = make_index(ref, *args)
 
        sys.stderr.write("Mapping sequences\n")
        reads = []
        scores = []

        extender = make_extender()
        for sequences in tqdm.tqdm(stream_fastq_pair(f1, f2)):
            r1, r2 = sequences
            a1 = extender.extend(index, r1)
            a2 = extender.extend(index, r2)
            reads.append(r1.id)
            reads.append(r2.id)
            #scores.append(score(index.fetch_ref(sam_ref[0].reference_start, sam_ref[0].reference_end), index.fetch_ref(a1.ref_start, a1.ref_end)))
            #scores.append(score(index.fetch_ref(sam_ref[1].reference_start, sam_ref[1].reference_end), index.fetch_ref(a2.ref_start, a2.ref_end)))
        reads = np.array(reads)
        scores = np.array(scores)
        #print("Mean edit distance between aligned and true sequence:", scores.mean())
        print("DONE")

from my_aligner import *

eval(FMIndexer, MyExtender)
#eval(HashIndexer, NaiveExtender, 19)