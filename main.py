# Import NeedlemanWunsch class and read_fasta function
from align import read_fasta, NeedlemanWunsch

def main():
    """
    This function should
    (1) Align all species to humans and print species in order of most similar to human BRD
    (2) Print all alignment scores between each species BRD2 and human BRD2
    """
    hs_seq, hs_header = read_fasta("./data/Homo_sapiens_BRD2.fa")
    gg_seq, gg_header = read_fasta("./data/Gallus_gallus_BRD2.fa")
    mm_seq, mm_header = read_fasta("./data/Mus_musculus_BRD2.fa")
    br_seq, br_header = read_fasta("./data/Balaeniceps_rex_BRD2.fa")
    tt_seq, tt_header = read_fasta("./data/tursiops_truncatus_BRD2.fa")

    # TODO Align all species to humans and print species in order of most similar to human BRD
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix
    alignment_dict = {}
    sequence_dict = {}

    NW_initialization = NeedlemanWunsch(sub_matrix_file = 'BLOSUM62.mat', gap_open = -10, gap_extend = -1)
    gallus_gallus_ = NW_initialization.align(seqA = hs_seq, seqB = gg_seq)
    gg_align_score = gallus_gallus_[0]
    alignment_dict[gg_align_score] = "Gallus Gallus"
    sequence_dict["Gallus Gallus"] = gallus_gallus_[2]
    

    mus_musculus_ = NW_initialization.align(seqA = hs_seq, seqB = mm_seq)
    mm_align_score = mus_musculus_[0]
    alignment_dict[mm_align_score] = "Mus Musculus"
    sequence_dict["Mus Musculus"] = mus_musculus_[2]

    balaeniceps_rex_ = NW_initialization.align(seqA = hs_seq, seqB = br_seq)
    br_align_score = balaeniceps_rex_[0]
    alignment_dict[br_align_score] = "Balaeniceps Rex"
    sequence_dict["Balaeniceps Rex"] = balaeniceps_rex_[2]

    tursiops_trancatus_ = NW_initialization.align(seqA = hs_seq, seqB = tt_seq)
    tt_align_score = tursiops_trancatus_[0]
    alignment_dict[tt_align_score] = "Tursiops Trancatus"
    sequence_dict["Tursiops Trancatus"] = tursiops_trancatus_[2]


    sorted_alignment_dict = sorted(alignment_dict.items(), reverse=True)    #ordered greatest score to least score

    for key, value in sorted_alignment_dict.items():
        return value, key   # returns the sequence(aligned to human) and species

        for k, v in sequence_dict.items():
            if value == k:
                return v  # returns alignment score

    # TODO print all of the alignment score between each species BRD2 and human BRD2
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix
    
    # See above!

if __name__ == "__main__":
    main()
