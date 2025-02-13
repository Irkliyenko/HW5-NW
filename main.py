# Import NeedlemanWunsch class and read_fasta function
from align import read_fasta, NeedlemanWunsch

def main():
    """
    This function should
    (1) Align all species to humans and print species in order of most similar to human BRD
    (2) Print all alignment scores between each species BRD2 and human BRD2
    """
    hs_seq, hs_header = read_fasta("./data/Homo_sapiens_BRD2.fa")

    # Create a dictionary where the key is the species name and the value is its sequence
    species_data = {
        "Gallus gallus": read_fasta("./data/Gallus_gallus_BRD2.fa")[0],
        "Mus musculus": read_fasta("./data/Mus_musculus_BRD2.fa")[0],
        "Balaeniceps rex": read_fasta("./data/Balaeniceps_rex_BRD2.fa")[0],
        "Tursiops truncatus": read_fasta("./data/tursiops_truncatus_BRD2.fa")[0]
    }

    # TODO Align all species to humans and print species in order of most similar to human BRD
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix
    sub_mat = './substitution_matrices/BLOSUM62.mat'

    # Initialize the Needleman-Wunsch class
    al = NeedlemanWunsch(sub_mat, gap_open = -10, gap_extend = -1)

    # Initialize a dictionary to store species and their alignment scores
    aln_score = {}

    # Align each species to human sequence
    for species, seq in species_data.items():
        score, aligned_hs, aligned_species = al.align(hs_seq, seq)
        aln_score[species] = score

    # Sort the dictionary by alignment score in descending order (most similar species first)
    sorted_dict = sorted(aln_score.items(), key=lambda x:x[1], reverse=True)

    # Print the species name and its alignment score
    for key, v in sorted_dict:
        print(f"Species: {key}\nAlignment score: {v}\n")


if __name__ == "__main__":
    main()
