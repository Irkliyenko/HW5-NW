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
    sub_mat = './substitution_matrices/BLOSUM62.mat'

    # Initialize the Needleman-Wunsch class
    nw = NeedlemanWunsch(sub_mat, gap_open=-10, gap_extend=-1)

    # Perform alignment
    score, seq1_aln, seq2_aln = nw.align(seq1, seq2)


    # Test if matrices are initialized
    assert nw._align_matrix is not None, "Alignment matrix wasn't initialized"
    assert nw._gapA_matrix is not None, "GapA matrix wasn't initialized"
    assert nw._gapB_matrix is not None, "GapB matrix wasn't initialized"

    # Test if matrices have the correct dimensions
    expected_shape = (len(seq1) + 1, len(seq2) + 1)
    assert nw._align_matrix.shape == expected_shape, f"Alignment matrix has incorrect dimensions: {nw._align_matrix.shape}"
    assert nw._gapA_matrix.shape == expected_shape, f"GapA matrix has incorrect dimensions: {nw._gapA_matrix.shape}"
    assert nw._gapB_matrix.shape == expected_shape, f"GapB matrix has incorrect dimensions: {nw._gapB_matrix.shape}"

    # Test if the first cell of the alignment matrix is correctly initialized
    assert nw._align_matrix[0, 0] == 0, "Top-left of the alignment matrix should be 0"

    # Test if the first row and column are initialized with correct gap penalties
    for i in range(1, len(seq1) + 1):
        expected_gap_score = nw.gap_open + (i) * nw.gap_extend
        assert nw._align_matrix[i, 0] == expected_gap_score, f"Incorrect gap penalty at row {i}, col 0"

    for j in range(1, len(seq2) + 1):
        expected_gap_score = nw.gap_open + (j) * nw.gap_extend
        assert nw._align_matrix[0, j] == expected_gap_score, f"Incorrect gap penalty at row 0, col {j}"

        assert nw._align_matrix[0, j] == expected_gap_score, f"Incorrect gap penalty at row 0, col {j}"


    

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
    sub_mat = './substitution_matrices/BLOSUM62.mat'

    # Initialize the Needleman-Wunsch class
    nw = NeedlemanWunsch(sub_mat, gap_open=-10, gap_extend=-1)

    # Perform alignment
    score, seq3_aln, seq4_aln = nw.align(seq3, seq4)
    print(score)
    print(seq3_aln)
    print(seq4_aln)

    # Test if the alignment output has the correct data types
    assert isinstance(score, int), "Score should be a float"
    assert isinstance(seq3_aln, str), "Aligned sequence 3 should be a string"
    assert isinstance(seq4_aln, str), "Aligned sequence 4 should be a string"

    # Test if the output sequences are not empty
    assert seq3_aln.strip() != "", "Aligned sequence 3 is empty"
    assert seq4_aln.strip() != "", "Aligned sequence 4 is empty"

    # Test if the sequences after backtracing have the same length and follow expected alignment rules
    assert len(seq3_aln) == len(seq4_aln), "Aligned sequences should have the same length"
    #assert all(c1 == '-' or c2 == '-' or c1 == c2 for c1, c2 in zip(seq3_aln, seq4_aln)), "Alignment does not follow expected behavior"

    # Test the score
    assert score == 17, "Alignment score isn't correct"
