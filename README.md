# Project 5
Needleman Wunsch Algorithm


# Project description
For this project, I completed the `NeedlemanWunsch.align` method, which implements the Needleman-Wunsch alignment algorithm.

First, the function initializes empty matrices to store matching scores, penalty scores, and a backtrace matrix that records the alignment path. If the maximum score at a given position equals the calculated match score, we move diagonally ("D"). If the maximum score corresponds to the gapA score, we move up ("U"). If it corresponds to the gapB score, we move left ("L").

Next, I completed the `NeedlemanWunsch._backtrace` method, which follows the backtrace matrix and returns the optimal alignment and alignment score.

I also wrote unit tests to verify the correctness of the `NeedlemanWunsch.align` and `NeedlemanWunsch._backtrace` methods.

The `test_nw_alignment` function verifies that the alignment matrices are correctly initialized and filled. It ensures that:
	▪️ The alignment matrix (`_align_matrix`) and gap matrices (`_gapA_matrix`, `_gapB_matrix`) are properly created.
	▪️ The matrices have the expected dimensions based on the input sequences.
	▪️ The first row and column are correctly initialized with gap penalties.
	▪️ The top-left cell of the alignment matrix is set to 0.

The `test_nw_backtrace` function ensures that the backtrace correctly reconstructs the optimal alignment. It checks that:
	▪️ The alignment score is correct and has the expected data type.
	▪️ The aligned sequences are non-empty and of equal length.
	▪️ The alignment score matches the expected value (17).
	▪️ The final aligned sequences match the expected output.
These tests help confirm that the Needleman-Wunsch algorithm correctly aligns sequences and accurately reconstructs the optimal alignment path.

