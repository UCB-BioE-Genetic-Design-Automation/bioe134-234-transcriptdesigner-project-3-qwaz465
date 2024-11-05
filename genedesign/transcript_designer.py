from genedesign.rbs_chooser import RBSChooser
from genedesign.models.transcript import Transcript
import csv
import random
from genedesign.checkers.codon_checker_window import CodonCheckerWindow
from genedesign.checkers.forbidden_sequence_checker import ForbiddenSequenceChecker
from genedesign.checkers.hairpin_checker import *
from genedesign.checkers.internal_promoter_checker import PromoterChecker
from genedesign.checkers.internal_rbs_checker import InternalRBSChecker

class TranscriptDesigner:
    """
    Reverse translates a protein sequence into a DNA sequence and chooses an RBS using the highest CAI codon for each amino acid.
    """

    def __init__(self):
        self.codonUsage = {}
        self.rbsChooser = None
        self.codonLists = {}
        self.codonChecker = None
        self.forbiddenSeqChecker = None
        self.PromoterChecker = None
        self.rbsChecker = None
        random.seed(1738)
    def initiate(self) -> None:
        """
        Initializes the codon map and the RBS chooser.
        """
        self.rbsChooser = RBSChooser()
        self.rbsChooser.initiate()
        # Open the file and read its contents
        with open('genedesign/data/codon_usage.txt', 'r') as file:
            reader = csv.reader(file, delimiter='\t')

            for row in reader:
                codon, amino_acid, relative_frequency, *_ = row
                # Store the relative frequency for each codon under its corresponding amino acid
                if amino_acid not in self.codonUsage:
                    self.codonUsage[amino_acid] = {}
        
                self.codonUsage[amino_acid][codon] = float(relative_frequency)
        self._generate_codon_lists()
        # boot up checkers
        self.codonChecker = CodonCheckerWindow()
        self.codonChecker.initiate()
        self.forbiddenSeqChecker = ForbiddenSequenceChecker()
        self.forbiddenSeqChecker.initiate()
        self.PromoterChecker = PromoterChecker()
        self.PromoterChecker.initiate()
        self.rbsChecker = InternalRBSChecker()
        self.rbsChecker.initiate()

    def run(self, peptide: str, ignores: set) -> Transcript:
        """
        Translates the peptide sequence to DNA and selects an RBS.
        
        Parameters:
            peptide (str): The protein sequence to translate.
            ignores (set): RBS options to ignore.
        
        Returns:
            Transcript: The transcript object with the selected RBS and translated codons.
        """
        # TODO add 5th checker  
        # TODO only check a few AAs back for preamble
        # put codons in here and then ''.join() when needed
        #should seed be here?
        codons = self.get_initial_cds(peptide)
        # TODO implement sliding window search with prev codons, curr window, and downstream stuff, already set up monster checker
        length = len(peptide)
        # TODO test get_initial_cds
        for window_start in range(3, length - 9 + 1, 3):
            # Get the current window and the next 6 amino acids (total of 9 amino acids)
            aa_window = peptide[window_start:window_start + 9]

            # Optimize codons for the current window using sliding window approach, considering upstream codons
            optimized_codons = self.sliding_window(codons, aa_window)
            codons.extend(optimized_codons)
            #TODO look at step 3 in gpt chat and add it here
            # take difference in list sizes for remaining_peptide then check indexs at the end
        remaining_peptide_index = len(codons)  # Remaining amino acids after full windows
        print(remaining_peptide_index) # this is how many have been done
        if remaining_peptide_index:
            # Use optimize_remaining for the last chunk of amino acids
            remaining_peptide = peptide[remaining_peptide_index:]
            remaining_codons = self.optimize_remaining(remaining_peptide)
            codons.extend(remaining_codons)
        print(f"starting pep length {len(peptide)} and end length {len(codons)}")
        # adds 6 to transcript somehow, hopefully through the rbs???
        stop_codon = self.guided_random('*')
        codons.extend(stop_codon)
        cds = ''.join(codons)

        # Choose an RBS
        selectedRBS = self.rbsChooser.run(cds, ignores)

        # Return the Transcript object
        return Transcript(selectedRBS, peptide, codons)
    
    def _generate_codon_lists(self):
        """
        Creates a 100-element list of codons for each amino acid where codons appear proportional to their relative frequency.
        """
        for amino_acid, codons in self.codonUsage.items():
            codon_list = []

            for codon, frequency in codons.items():
                # Calculate how many times the codon should appear in a 100-element list
                count = int(round(frequency * 100))
                codon_list.extend([codon] * count)

            # If the codon list has more than 100 elements (due to rounding), trim it
            if len(codon_list) > 100:
                codon_list = codon_list[:100]

            # If the codon list has fewer than 100 elements, pad it with the most frequent codon
            elif len(codon_list) < 100:
                most_frequent_codon = max(codons, key=codons.get)  # Find the most frequent codon
                codon_list.extend([most_frequent_codon] * (100 - len(codon_list)))

            # Store the codon list for the amino acid
            self.codonLists[amino_acid] = codon_list
    
    def guided_random(self, aa : str) -> str:
        # setting seed here, then whats the point of sampling more than once since im gonna get the same codon each time i run this???
        codon = random.choice(self.codonLists[aa])
        return codon
    
    def get_initial_cds(self, peptide : str):
        #TODO should have no error, return best one
        # currently picks the one that passes the most checkers, see if i can potentially optimize using checker values?
        best_codons = ''
        most_checkers = -1
        for y in range(100):
            print(f"initial cds attempt {y}")
            codons = [self.guided_random(peptide[x]) for x in range(3)]
            print(codons)
            checker_passes = self.monster_checker(codons)
            print(f"checkers passed: {checker_passes}")
            # print(good_cds)
            #TODO change this if adding checkers
            if checker_passes == 5:
                return codons
            if checker_passes > most_checkers:
                most_checkers = checker_passes
                best_codons = codons
            return best_codons

            
    def sliding_window(self, upstream_codons, aa_window):
        """
        Optimizes codons for the first 3 amino acids in the given window, considering upstream codons and the full window (9 amino acids).
        
        Parameters:
            upstream_codons (list): List of previously chosen codons (upstream context).
            aa_window (str): A 9-amino-acid window (the first 3 are the current window, the remaining 6 are context).
        
        Returns:
            list: A list of optimized codons for the first 3 amino acids.
        """
        best_codons = ''
        most_checkers = -1
        window_size = 3
        downstream_aa_count = 6

        # Ensure the upstream codons are at least 6 long, pad if necessary
        upstream_window_start = max(0, len(upstream_codons) - 6)
        codons_before = upstream_codons[upstream_window_start:]

        # Optimize codons for the current 9-amino-acid window
        for attempt in range(100):  # Try 100 times to find the best codons
            print(f"Sliding window attempt {attempt}")
            current_codons = codons_before.copy()  # Include the 6 upstream codons for context
            for x in current_codons:
                print(f"previous codon : {x}")

            # Optimize codons for the 9 amino acids in the window (first 3 are the current window)
            window_codons = []
            for aa in aa_window[:window_size]:  # Optimize for the first 3 amino acids
                codon = self.guided_random(aa)
                print(f"Window codon: {codon}")
                window_codons.append(codon)

            # Generate codons for the downstream 6 amino acids
            downstream_codons = []
            for aa in aa_window[window_size:window_size + downstream_aa_count]:
                codon = self.guided_random(aa)
                print(f"Downstream codon: {codon}")
                downstream_codons.append(codon)

            # Concatenate codons before, window codons, and downstream codons
            full_sequence = current_codons + window_codons + downstream_codons

            # Check the entire codon sequence (6 upstream + window + 6 downstream) using monster_checker
            checker_passes = self.monster_checker(full_sequence)
            print(f"Checkers passed: {checker_passes}")

            # If all checks pass, return the first 3 codons of the current window
            if checker_passes == 5:  # Maximum score assumed to be 5
                return window_codons

            # Keep track of the best codons if no perfect match
            if checker_passes > most_checkers:
                most_checkers = checker_passes
                best_codons = window_codons

        # Return the best codons for the first 3 amino acids if no perfect match was found
        return best_codons
    
    def optimize_remaining(self, remaining_peptide):
        """
        Optimizes the codons for the remaining amino acids at the end of the peptide sequence.

        Parameters:
            prev_codons (list): Codons translated before this remaining sequence.
            remaining_peptide (str): The final sequence of amino acids (less than 9).

        Returns:
            list: Optimized codons for the remaining amino acids.
        """
        codons = []
        best_codons = ''
        most_checkers = -1

        # Handle the remaining amino acids   
        for x in range(100):
            print(f"remaining attempt {x} with remaining length {len(remaining_peptide)}")
            codons = [self.guided_random(aa) for aa in remaining_peptide]
            print(codons)

            # Run the monster_checker to ensure the remaining sequence passes the checks
            checker_passes = self.monster_checker(codons)
            print(f"checkers passed: {checker_passes}")
            # If all checks pass, return the first 3 codons of the current window
            if checker_passes == 5:  # Maximum score assumed to be 5
                return codons  # Return first 3 codons of the current window
                # Keep track of the best codons if no perfect match
            if checker_passes > most_checkers:
                most_checkers = checker_passes
                best_codons = codons

        # Return the best codons for the first 3 amino acids if no perfect match was found
        return best_codons

    def monster_checker(self, cds):
        #TODO also return specific values 
        cds_str = ''.join(cds)
        # bool is 0th index in returned tuple
        codon_checker_result = self.codonChecker.run(cds)
        codon_checker_pass =codon_checker_result[0]
        cai = codon_checker_result[3] 
        forbidden_seq_checker_result = self.forbiddenSeqChecker.run(cds_str)
        forbidden_seq_pass = forbidden_seq_checker_result[0]
        hairpin_checker_result = hairpin_checker(cds_str)
        hairpin_pass = hairpin_checker_result[0]
        promoter_checker_pass = self.PromoterChecker.run(cds_str)[0]
        rbs_checker_pass = self.rbsChecker.run(cds_str)[0]
        passes = [codon_checker_pass, forbidden_seq_pass, hairpin_pass, promoter_checker_pass, rbs_checker_pass]
        checker_names = ['cai checker', 'forbidden seq checker', 'hairpin checker', 'promoter checker', 'rbs checker']
        for x in range(5):
            if passes[x] == False:
                print(checker_names[x])
        return codon_checker_pass + forbidden_seq_pass + hairpin_pass + promoter_checker_pass + rbs_checker_pass
if __name__ == "__main__":
    # Example usage of TranscriptDesigner
    peptide = "MYPFIRTARMTV"
    
    designer = TranscriptDesigner()
    designer.initiate()

    ignores = set()
    transcript = designer.run(peptide, ignores)
    
    # Print out the transcript information
    print(transcript)
