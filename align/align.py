# Importing Dependencies
import numpy as np
from typing import Tuple

# Defining class for Needleman-Wunsch Algorithm for Global pairwise alignment
class NeedlemanWunsch:
    """ Class for NeedlemanWunsch Alignment

    Parameters:
        sub_matrix_file: str
            Path/filename of substitution matrix
        gap_open: float
            Gap opening penalty
        gap_extend: float
            Gap extension penalty

    Attributes:
        seqA_align: str
            seqA alignment
        seqB_align: str
            seqB alignment
        alignment_score: float
            Score of alignment from algorithm
        gap_open: float
            Gap opening penalty
        gap_extend: float
            Gap extension penalty
    """
    def __init__(self, sub_matrix_file: str, gap_open: float, gap_extend: float):
        # Init alignment and gap matrices
        self._align_matrix = None
        self._gapA_matrix = None
        self._gapB_matrix = None

        # Init matrices for backtrace procedure
        self._back = None
        self._back_A = None
        self._back_B = None

        # Init alignment_score
        self.alignment_score = 0

        # Init empty alignment attributes
        self.seqA_align = ""
        self.seqB_align = ""

        # Init empty sequences
        self._seqA = ""
        self._seqB = ""

        # Setting gap open and gap extension penalties
        self.gap_open = gap_open
        assert gap_open < 0, "Gap opening penalty must be negative."
        self.gap_extend = gap_extend
        assert gap_extend < 0, "Gap extension penalty must be negative."

        # Generating substitution matrix
        self.sub_dict = self._read_sub_matrix(sub_matrix_file) # substitution dictionary

    def _read_sub_matrix(self, sub_matrix_file):
        """
        DO NOT MODIFY THIS METHOD! IT IS ALREADY COMPLETE!

        This function reads in a scoring matrix from any matrix like file.
        Where there is a line of the residues followed by substitution matrix.
        This file also saves the alphabet list attribute.

        Parameters:
            sub_matrix_file: str
                Name (and associated path if not in current working directory)
                of the matrix file that contains the scoring matrix.

        Returns:
            dict_sub: dict
                Substitution matrix dictionary with tuple of the two residues as
                the key and score as value e.g. {('A', 'A'): 4} or {('A', 'D'): -8}
        """
        with open(sub_matrix_file, 'r') as f:
            dict_sub = {}  # Dictionary for storing scores from sub matrix
            residue_list = []  # For storing residue list
            start = False  # trigger for reading in score values
            res_2 = 0  # used for generating substitution matrix
            # reading file line by line
            for line_num, line in enumerate(f):
                # Reading in residue list
                if '#' not in line.strip() and start is False:
                    residue_list = [k for k in line.strip().upper().split(' ') if k != '']
                    start = True
                # Generating substitution scoring dictionary
                elif start is True and res_2 < len(residue_list):
                    line = [k for k in line.strip().split(' ') if k != '']
                    # reading in line by line to create substitution dictionary
                    assert len(residue_list) == len(line), "Score line should be same length as residue list"
                    for res_1 in range(len(line)):
                        dict_sub[(residue_list[res_1], residue_list[res_2])] = float(line[res_1])
                    res_2 += 1
                elif start is True and res_2 == len(residue_list):
                    break
        return dict_sub

    def match(self, residueA: str, residueB: str) -> int:
        """
        Creating a function to determine the match score between two amino acid residues

        Parameters:
            residueA: str
                string of a single residue to be matched in the substitution matrix dictionary

            residueB: str
                string of the other single residue to be matched in the substitution matrix dictionary

        Returns:
            score: int
                the match/mismatch score between residueA and residueB
        """
        dict_sub = self.sub_dict

        for key, value in dict_sub.items():
            if key == (residueA, residueB):
                score = value

        return int(score)

    def align(self, seqA: str, seqB: str) -> Tuple[float, str, str]:
        """
        # TODO: Fill in the Needleman-Wunsch Algorithm below
        to perform global sequence alignment of seqA and seqB
        and return a tuple with the following format
        (alignment score, seqA alignment, seqB alignment)

        This function takes 2 sequences and determines the best alignment (with or without
        gaps in the sequence). It uses the match() function to determine match score from a 
        substitution scoring matrix.

        Parameters:
            seqA: str
                string of an amino acid residue sequence
            seqB: str
                string of another amino acid residue sequence to be aligned with seqA

        Returns:
            alignment_score, seqA_align, seqB_align : tuple
                tuple of the alignment score, with its respective sequence A alignment and sequence B alignment
        """
        # Initialize 6 matrix private attributes for use in alignment
        # create matrices for alignment scores and gaps
        self._align_matrix = np.ones((len(seqA) + 1, len(seqB) + 1)) * -np.inf
        self._gapA_matrix = np.ones((len(seqA) + 1, len(seqB) + 1)) * -np.inf
        self._gapB_matrix = np.ones((len(seqA) + 1, len(seqB) + 1)) * -np.inf

        # create matrices for pointers used in backtrace procedure
        self._back = np.ones((len(seqA) + 1, len(seqB) + 1)) * -np.inf
        self._back_A = np.ones((len(seqA) + 1, len(seqB) + 1)) * -np.inf
        self._back_B = np.ones((len(seqA) + 1, len(seqB) + 1)) * -np.inf

        # Resetting alignment in case method is called more than once
        self.seqA_align = ""
        self.seqB_align = ""

        # Resetting alignment score in case method is called more than once
        self.alignment_score = 0

        # Initializing sequences for use in backtrace method
        self._seqA = seqA
        self._seqB = seqB

        # TODO Implement the global sequence alignment here

        # initialize first row and first column of each matrix
        # initialize backtracing matrices as well (O = up, 1 = diagonal, 2 = left)

        self._align_matrix[0][0] = 0
        self._gapB_matrix[:, 0] = self.gap_open + self.gap_extend * np.array(range(len(self._seqA) + 1))
        self._gapA_matrix[0,:] = self.gap_open + self.gap_extend * np.array(range(len(self._seqB) + 1))

        for i in range(1, len(seqA)+1):
            self._back_B[i][0] = 0

        for j in range(1, len(seqB)+1):
            self._back_A[0][j] = 2

        # filling in the matrices row by row
        M = self._align_matrix
        Ix = self._gapB_matrix       #along the x axis, is sequence B
        Iy = self._gapA_matrix       #along the y axis, is sequence A

        for i in range(1, len(seqA)+1):
            for j in range(1, len(seqB)+1):

                # alignment matrix (match or mismatch)
                M_M = self._align_matrix[i-1][j-1] + self.match(seqA[i-1], seqB[j-1])
                M_Ix = self._gapB_matrix[i-1][j-1] + self.match(seqA[i-1], seqB[j-1])
                M_Iy = self._gapA_matrix[i-1][j-1] + self.match(seqA[i-1], seqB[j-1])

                # gapB matrix (gaps in seqB)
                Ix_M = self._align_matrix[i-1][j] + self.gap_extend + self.gap_open
                Ix_Ix = self._gapB_matrix[i-1][j] + self.gap_extend
                Ix_Iy = self._gapA_matrix[i-1][j] + self.gap_extend + self.gap_open

                # gapA matrix (gaps in seqA)
                Iy_M = self._align_matrix[i][j-1] + self.gap_extend + self.gap_open
                Iy_Ix = self._gapB_matrix[i][j-1] + self.gap_extend + self.gap_open
                Iy_Iy = self._gapA_matrix[i][j-1] + self.gap_extend

                # appending maximum value to matrices
                # Adding code to backtracing matrices to tell you where to move next (highroad alignment)
                # 0: if Ix was the highest (gap in seqB), move [i-1][j]
                # 1: if M was the highest (match/mismatch), move [i-1][j-1]
                # 2: if Iy was the highest (gap is seqA), move [i][j-1])

                max_value_M = max(M_Ix, M_M, M_Iy)
                self._align_matrix[i][j] = max_value_M
                self._back[i][j] = [M_Ix, M_M, M_Iy].index(max_value_M)

                max_value_Ix = max(Ix_Ix, Ix_M, Ix_Iy)
                self._gapB_matrix[i][j] = max_value_Ix
                self._back_B[i][j] = [Ix_Ix, Ix_M, Ix_Iy].index(max_value_Ix) 


                max_value_Iy = max(Iy_Ix, Iy_M, Iy_Iy)
                self._gapA_matrix [i][j] = max_value_Iy
                self._back_A[i][j] = [Iy_Ix, Iy_M, Iy_Iy].index(max_value_Iy) 
          


        return self._backtrace()

    def _backtrace(self) -> Tuple[float, str, str]:
        """
        # TODO Implement the traceback procedure method below
        based on the heuristic you implement in the align method.
        The traceback method should return a tuple of the alignment
        score, the seqA alignment and the seqB alignment respectively.
        """
        # Implementing this method based upon the heuristic chosen in the align method above (highroad alignment).
        M = self._align_matrix
        Ix = self._gapB_matrix
        Iy = self._gapA_matrix

        last_col = len(self._seqB)
        last_row = len(self._seqA)

        self.alignment_score = max(Ix[last_row][last_col], M[last_row][last_col], Iy[last_row][last_col])

        # creating objects for the for loop
        _row = len(self._seqA)
        _column = len(self._seqB)

        seqA_align = ""
        seqB_align = ""

        # location of the maximum value (in which matrix)
        curr_location = [Ix[last_row][last_col], M[last_row][last_col], Iy[last_row][last_col]].index(self.alignment_score)

        while _row > 0 or _column > 0:

            pointer = curr_location 

            # gapB matrix
            if pointer == 0:
                seqA_align += self._seqA[_row-1]
                seqB_align += '-'
                curr_location = self._back_B[_row][_column]
                _row = _row -1
                
            #alignment matrix
            elif pointer == 1:
                seqA_align += self._seqA[_row-1]
                seqB_align += self._seqB[_column-1]
                curr_location = self._back[_row][_column]
                _row = _row - 1
                _column = _column - 1

            # gapA matrix
            elif pointer == 2:
                seqA_align += '-'
                seqB_align += self._seqB[_column-1]
                curr_location = self._back_A[_row][_column]
                _column = _column - 1


        final_seqA_align = seqA_align[::-1]
        final_seqB_align = seqB_align[::-1]
        
        return self.alignment_score, final_seqA_align, final_seqB_align




def read_fasta(fasta_file: str) -> Tuple[str, str]:
    """
    DO NOT MODIFY THIS FUNCTION! IT IS ALREADY COMPLETE!

    This function reads in a FASTA file and returns the associated
    string of characters (residues or nucleotides) and the header.
    This function assumes a single protein or nucleotide sequence
    per fasta file and will only read in the first sequence in the
    file if multiple are provided.

    Parameters:
        fasta_file: str
            name (and associated path if not in current working directory)
            of the Fasta file.

    Returns:
        seq: str
            String of characters from FASTA file
        header: str
            Fasta header
    """
    assert fasta_file.endswith(".fa"), "Fasta file must be a fasta file with the suffix .fa"
    with open(fasta_file) as f:
        seq = ""  # initializing sequence
        first_header = True
        for line in f:
            is_header = line.strip().startswith(">")
            # Reading in the first header
            if is_header and first_header:
                header = line.strip()  # reading in fasta header
                first_header = False
            # Reading in the sequence line by line
            elif not is_header:
                seq += line.strip().upper()  # generating full sequence
            # Breaking if more than one header is provided in the fasta file
            elif is_header and not first_header:
                break
    return seq, header
