# Importing Dependencies
import pytest
from align import NeedlemanWunsch, read_fasta
import numpy as np

def test_nw_alignment():
    """
    TODO: Write your unit test for NW alignment
    using test_seq1.fa and test_seq2.fa by
    asserting that you have correctly filled out
    the your 3 alignment matrices.
    Use the BLOSUM62 matrix and a gap open penalty
    of -10 and a gap extension penalty of -1.
    """
    seq1, _ = read_fasta("./data/test_seq1.fa")
    seq2, _ = read_fasta("./data/test_seq2.fa")
    
    N_W = NeedlemanWunsch(sub_matrix_file = "./substitution_matrices/BLOSUM62.mat", gap_open = -10, gap_extend = -1)
    N_W_align = N_W.align(seqA = seq1, seqB = seq2)

    # test alignment accuracy here: check each matrix has same dimensions and there are no null values?

    assert N_W._back.shape == N_W._back_A.shape

    for i in range(len(N_W._back.shape[0])):
        for j in range(len(N_W._back.shape[1])):
            assert N_W._align_matrix[i][j] != None


def test_nw_backtrace():
    """
    TODO: Write your unit test for NW backtracing
    using test_seq3.fa and test_seq4.fa by
    asserting that the backtrace is correct.
    Use the BLOSUM62 matrix. Use a gap open
    penalty of -10 and a gap extension penalty of -1.
    """
    seq3, _ = read_fasta("./data/test_seq3.fa")
    seq4, _ = read_fasta("./data/test_seq4.fa")

    N_W = NeedlemanWunsch(sub_matrix_file = "./substitution_matrices/BLOSUM62.mat", gap_open = -10, gap_extend = -1)
    N_W_align = N_W.align(seqA = seq3, seqB = seq4)

    # test backtracing accuracy here: check that sequence score (manually scoring) = alignment score

    alignment_score = N_W_align[0]
    sequence_A = N_W_align[1]
    sequence_B = N_W_align[2]

    sequence_score = 0

    for i in range(len(sequence_A)):
        if sequence_A[i] != '-':
            if sequence_B[i] != '-':
                for key, value in N_W.sub_dict.items():
                    if key == (sequence_A[i], sequence_B[i]):
                        sequence_score += value
            elif sequence_B[i] == '-':
                if i == 0:
                    sequence_score += -11
                elif i > 0:
                    if sequence_B[i-1] != '-':
                        sequence_score += -11
                    elif sequence_B[i-1] == '-':
                        sequence_score += -1

        elif sequence_A[i] == '-':
            if i == 0:
                sequence_score += -11
            elif i > 0 :
                if sequence_A[i-1] != '-':
                    sequence_score += -11
                elif sequence_A[i-1] == '-':
                    sequence_score += -1

    assert alignment_score == sequence_score





