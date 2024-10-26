from genedesign.rbs_chooser import RBSChooser
from genedesign.models.transcript import Transcript
import csv
import random

class TranscriptDesigner:
    """
    Reverse translates a protein sequence into a DNA sequence and chooses an RBS using the highest CAI codon for each amino acid.
    """

    def __init__(self):
        self.codonUsage = {}
        self.rbsChooser = None
        self.codonLists = {}

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

    def run(self, peptide: str, ignores: set) -> Transcript:
        #TODO implement monte carlo selection for first group, some of this may need to be done in initiate?
        """
        Translates the peptide sequence to DNA and selects an RBS.
        
        Parameters:
            peptide (str): The protein sequence to translate.
            ignores (set): RBS options to ignore.
        
        Returns:
            Transcript: The transcript object with the selected RBS and translated codons.
        """


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

if __name__ == "__main__":
    # Example usage of TranscriptDesigner
    peptide = "MYPFIRTARMTV"
    
    designer = TranscriptDesigner()
    designer.initiate()

    ignores = set()
    transcript = designer.run(peptide, ignores)
    
    # Print out the transcript information
    print(transcript)
