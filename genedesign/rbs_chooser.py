from genedesign.models.rbs_option import RBSOption
from genedesign.rbs_chooser_data_cleaning import *
from genedesign.seq_utils.Translate import *
from genedesign.seq_utils.hairpin_counter import *
from genedesign.seq_utils.calc_edit_distance import *

class RBSChooser:
    """
    A simple RBS selection algorithm that chooses an RBS from a list of options, excluding any RBS in the ignore set.
    """

    def __init__(self):
        self.rbsOptions = []
        self.translator = Translate()
        self.translator.initiate()

    def initiate(self) -> None:
        """
        Initializes the internal state of the RBSChooser, including setting up the list of RBS options.

        This method loads RBSOptions into the `rbs_options_list` attribute by extracting data from
        the `genes_info` and `top_5_percent_list`. During this process, coding sequences (CDS) from
        gene data are translated into protein sequences to extract the first six amino acids for
        each RBS option.

        If an invalid CDS sequence is encountered during translation (e.g., non-nucleotide
        characters or improper length), a ValueError will be raised, indicating which locus tag
        caused the issue.

        Attributes initialized:
            rbs_options_list (list): A list of available RBS options, populated during initiation.

        Raises:
            ValueError: If an invalid CDS sequence is encountered during translation,
                        indicating the problematic locus tag.

        Authors:
        Jake Barer
        J. Christopher Anderson

        Version 1.0

        Last Modified 9/28/24
        """
        # data extraction
        genbank_file = "genedesign/data/sequence.gb"  # Now using the correct file name
        genes_info = extract_genes_info(genbank_file)
        file_path = 'genedesign/data/511145-WHOLE_ORGANISM-integrated.txt'
        top_5_percent_list = proteomics_prune(file_path)
        translator = self.translator
        # Loop through top_5_percent_list and genes_info to create RBSOption objects
        for locus_tag, abundance in top_5_percent_list:
            if locus_tag in genes_info:
                gene_info = genes_info[locus_tag]

                utr = gene_info['UTR']
                cds = gene_info['CDS']
                gene_name = gene_info['gene']

                # Translate the CDS to get the protein sequence
                try:
                    protein_sequence = translator.run(str(cds))  # Ensure the CDS is passed as a string
                    # Get the first six amino acids from the protein sequence
                    first_six_aas = protein_sequence[:6]

                    # Create an RBSOption instance
                    rbs_option = RBSOption(
                        utr=str(utr),  # Ensure the UTR is passed as a string
                        cds=str(cds),  # Ensure the CDS is passed as a string
                        gene_name=gene_name,
                        first_six_aas=first_six_aas
                    )

                    # Add the RBSOption instance to the class variable list
                    self.rbsOptions.append(rbs_option)

                except ValueError as e:
                    print(f"Error translating CDS for {locus_tag}: {e}")

    def run(self, cds: str, ignores: set) -> RBSOption:
        """
        Selects the best RBS option for a given coding sequence (CDS) using a combination
        of two metrics: the edit distance between the first six amino acids and the
        number of potential hairpin structures in the UTR + CDS region.

        Parameters:
            cds (str): The coding sequence (CDS) for which the best RBS is to be selected.
                       The first six amino acids of this CDS are compared with those of each RBS option.
            ignores (Set[RBSOption]): A set of RBSOption objects to ignore when searching for
                                      the best RBS option.

        Returns:
            RBSOption: The RBSOption that provides the best match to the input CDS based on
                       the lowest edit distance, and secondarily, the fewest potential hairpin structures.

        Raises:
            ValueError: If the input CDS is invalid (e.g., contains non-nucleotide characters),
                        if the CDS length is not a multiple of 3, if no valid RBS options remain
                        after filtering, or if the CDS sequence is empty.

        Authors:
            Jake Barer
            J. Christopher Anderson

        Version 1.0

        Last Modified 9/28/24
        """
        # Validate that the CDS is not empty
        if not cds:
            raise ValueError("CDS sequence cannot be empty.")

        # Validate that the CDS contains only valid nucleotide characters (A, T, C, G)
        if not all(base in 'ATCG' for base in cds):
            raise ValueError("Invalid CDS sequence: contains non-nucleotide characters. Only A, T, C, G are allowed.")

        # Validate that the CDS length is a multiple of 3 for proper translation
        if len(cds) % 3 != 0:
            raise ValueError("CDS sequence length must be a multiple of 3 for valid translation.")

        # Check if there are any RBS options to work with
        if not self.rbsOptions:
            raise ValueError("No RBS options are available to choose from.")

        # Filter out ignored RBS options
        valid_rbs_options = [rbs_option for rbs_option in self.rbsOptions if rbs_option not in ignores]
        if not valid_rbs_options:
            raise ValueError("No valid RBS options remain after applying the ignore filter.")

        # Initialize the Translate object and translate the first 6 amino acids of the input CDS
        translator = self.translator
        input_first_six_aas = translator.run(cds[:18])  # 18 bases for 6 amino acids

        best_rbs = None
        best_score = float('inf')  # Start with a very high score

        for rbs_option in valid_rbs_options:
            # Calculate the edit distance between input CDS first six AAs and RBS option first six AAs
            edit_distance = calculate_edit_distance(input_first_six_aas, rbs_option.first_six_aas)

            # Calculate the hairpin count for the UTR + CDS sequence of the RBS option
            hairpin_count = hairpin_counter(rbs_option.utr + rbs_option.cds)[0]

            # Combine the two scores; prioritize the edit distance, then hairpin count
            # doing it this way because if theres a high edit distance then there will be no binding so nothing will happen in the first place
            # however if theres a lot of secondary structure, translation can still occur but much worse, so I am placing much higher priority on
            # edit distance since that more or less says if binding occurs in the first place

            score = (edit_distance * 1000) + hairpin_count  # Giving much higher weight to edit distance

            # If this score is better, update the best RBS
            if score < best_score:
                best_score = score
                best_rbs = rbs_option

        return best_rbs
if __name__ == "__main__":
    # Example usage of RBSChooser
    cds = "ATGGTAAGAAAACAGTTGCAGAGAGTTGAATT..."

    # Initialize the chooser
    chooser = RBSChooser()
    chooser.initiate()

    # Choose RBS with no ignores
    ignores = set()
    selected1 = chooser.run(cds, ignores)
    
    # Add the first selection to the ignore list
    ignores.add(selected1)
    
    # Choose another RBS option after ignoring the first
    selected2 = chooser.run(cds, ignores)

    # Print the selected RBS options
    print("Selected1:", selected1)
    print("Selected2:", selected2)
