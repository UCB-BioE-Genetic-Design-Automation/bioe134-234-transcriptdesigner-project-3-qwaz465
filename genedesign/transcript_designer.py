from genedesign.rbs_chooser import RBSChooser
from genedesign.models.transcript import Transcript
import csv
import random
from genedesign.checkers.codon_checker import CodonChecker
from genedesign.checkers.forbidden_sequence_checker import ForbiddenSequenceChecker
from genedesign.checkers.hairpin_checker import *
from genedesign.checkers.internal_promoter_checker import PromoterChecker

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
        self.codonChecker = CodonChecker()
        self.codonChecker.initiate()
        self.forbiddenSeqChecker = ForbiddenSequenceChecker()
        self.forbiddenSeqChecker.initiate()
        self.PromoterChecker = PromoterChecker()
        self.PromoterChecker.initiate()

    def run(self, peptide: str, ignores: set) -> Transcript:
        """
        Translates the peptide sequence to DNA and selects an RBS.
        
        Parameters:
            peptide (str): The protein sequence to translate.
            ignores (set): RBS options to ignore.
        
        Returns:
            Transcript: The transcript object with the selected RBS and translated codons.
        """
        # put codons in here and then ''.join() when needed
        codons = self.get_initial_cds(peptide)
        # TODO implement sliding window search with prev codons, curr window, and downstream stuff, already set up monster checker
        # TODO test get_initial_cds


        # Choose an RBS
        # selectedRBS = self.rbsChooser.run(cds, ignores)

        # Return the Transcript object
        # return Transcript(selectedRBS, peptide, codons)
    
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
        random.seed(1738)
        codon = random.choice(self.codonLists[aa])
        return codon
    
    def get_initial_cds(self, peptide : str):
        #TODO should have no error, return best one
        good_cds = False
        for x in range(100):
            codons = []
            for x in range(2):
                aa = peptide[x]
                codon = self.guided_random(aa)
                codons.extend(codon)
            good_cds = self.monster_checker(codons)
            # print(good_cds)
            if good_cds:
                return codons
        raise RecursionError('no good cds present')
            
    
    def monster_checker(self, cds) -> bool:
        #TODO also return specific values 
        cds_str = ''.join(cds)
        # bool is 0th index in returned tuple
        codon_checker_result = self.codonChecker.run(cds)
        codon_checker_pass =codon_checker_result[0]
        cai = codon_checker_result[3] 
        forbidden_seq_checker_result = self.forbiddenSeqChecker.run(cds_str)
        forbidden_seq_result = forbidden_seq_checker_result[0]
        hairpin_checker_result = hairpin_checker(cds_str)[0]
        promoter_checker_result = self.PromoterChecker.run(cds_str)[0]
        
        if codon_checker_result and forbidden_seq_checker_result and hairpin_checker_result and promoter_checker_result:
            return True
        print(codon_checker_result)
        print(forbidden_seq_checker_result)
        print(hairpin_checker_result)
        print(promoter_checker_result)
        return False
if __name__ == "__main__":
    # Example usage of TranscriptDesigner
    peptide = "MYPFIRTARMTV"
    
    designer = TranscriptDesigner()
    designer.initiate()

    ignores = set()
    transcript = designer.run(peptide, ignores)
    
    # Print out the transcript information
    print(transcript)
