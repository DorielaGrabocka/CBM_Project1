import sys
import pickle
import itertools
from abc import ABC, abstractmethod
from typing import Generator, Optional

import pysam # type: ignore
from Bio import SeqIO, SeqRecord # type: ignore

## Storage classes
class Seed:
    ''' A storage class for seeds. '''
    query_id = ""
    query_seq: Optional[str] = None
    orientation: Optional[bool] = None
    query_start: Optional[int] = None
    query_end: Optional[int] = None
    ref_start: Optional[int] = None
    ref_end: Optional[int] = None

    def __init__(self,
                 query_id: str,
                 query_seq: Optional[str] = None,
                 orientation: Optional[bool] = None,
                 query_start: Optional[int] = None,
                 query_end: Optional[int] = None,
                 ref_start: Optional[int] = None,
                 ref_end: Optional[int] = None) -> None:
        self.query_id = query_id
        self.query_seq = query_seq
        self.orientation = orientation
        self.query_start = query_start
        self.query_end = query_end
        self.ref_start = ref_start
        self.ref_end = ref_end

    def __str__(self) -> str:
        return "{} {} {} [{},{}) [{},{})".format(
            self.query_id,
            self.query_seq if self.query_seq else "*",
            "*" if self.orientation is None else ("-" if self.orientation else "+"),
            "*" if self.query_start is None else self.query_start,
            "*" if self.query_end is None else self.query_end,
            "*" if self.ref_start is None else self.ref_start,
            "*" if self.ref_end is None else self.ref_end
        )

class Alignment(Seed):
    ''' A storage class for alignments. It additionally stores an alignment score
        and a CIGAR to indicate the edit path to the target sequence. '''
    score: Optional[int] = None
    cigar: Optional[str] = None

    def __init__(self,
                 query_id: str,
                 query_seq: Optional[str] = None,
                 orientation: Optional[bool] = None,
                 query_start: Optional[int] = None,
                 query_end: Optional[int] = None,
                 ref_start: Optional[int] = None,
                 ref_end: Optional[int] = None,
                 score: Optional[int] = None,
                 cigar: Optional[str] = None) -> None:
        super().__init__(query_id, query_seq, orientation, query_start, query_end, ref_start, ref_end)
        self.score = score
        self.cigar = cigar

    def __str__(self) -> str:
        return "{} {} {}".format(super().__str__(),
                                 "*" if self.score is None else self.score,
                                 "*" if self.cigar is None else self.cigar)

## Abstract parent classes
class ReferenceIndex(ABC):
    ''' A class for indexing the reference genome to allow for fast searches. '''
    ref_file: Optional[pysam.libcfaidx.FastaFile] = None

    def __init__(self, ref_file: str) -> None:
        try:
            self.load(ref_file + ".idx")
        except:
            sys.stderr.write(f'Failed to find index {ref_file}.idx, generating...\n')
            self.ref_file = pysam.FastaFile(ref_file)
            self._generate_index()

    @abstractmethod
    def _generate_index(self) -> None:
        ''' Construct the index '''
        pass

    def fetch_ref(self, begin: int, end: int) -> str:
        ''' Return the substring from the reference in the interval [begin, end) '''
        assert(not self.ref_file is None)
        return self.ref_file.fetch(reference=self.ref_file.references[0], start=begin, end=end)

    @abstractmethod
    def serialize(self, out_file: str) -> None:
        ''' Write the index to a file '''
        pass

    @abstractmethod
    def load(self, in_file: str) -> None:
        ''' Load the index from a file '''
        pass

    @abstractmethod
    def get_seeds(self, query_record: SeqRecord) -> list[Seed]:
        '''Given a query sequence, return all matching seeds in the index'''
        pass

class Extender(ABC):
    @abstractmethod
    def extend(self, index: ReferenceIndex, query_record: SeqRecord) -> Optional[Alignment]:
        ''' Return the best scoring alignment. This method should extract seeds using the index '''
        pass


## Helpers for reading query FASTQ files
def get_fastq_file_iterator(query_file: str) -> SeqIO.QualityIO.FastqPhredIterator:
    ''' Return an iterator for the sequences from a query FASTQ file '''
    return SeqIO.parse(query_file, "fastq")

def get_fasta_file_iterator(query_file: str):
    ''' Return an iterator for the sequences from a query FASTA file '''
    return SeqIO.parse(query_file, "fasta")

def stream_fastq_pair(fwd_file: str, rev_file: str):
    ''' Return an iterator for the forward and reverse reads in a paired read file '''
    return zip(get_fastq_file_iterator(fwd_file),
               get_fastq_file_iterator(rev_file))

def stream_fastq_ref_pair(fwd_file: str, rev_file: str, ref_file: str):
    ''' Return an iterator for the forward and reverse reads in a paired read file,
        and their corresponding reference sequences. '''
    return zip(get_fastq_file_iterator(fwd_file),
               get_fastq_file_iterator(rev_file),
               itertools.zip_longest(*[get_fasta_file_iterator(ref_file)]*2))

