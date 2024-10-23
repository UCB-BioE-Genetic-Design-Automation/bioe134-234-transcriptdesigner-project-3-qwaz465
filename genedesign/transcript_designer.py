from genedesign.rbs_chooser import RBSChooser
from genedesign.models.transcript import Transcript
import csv

class TranscriptDesigner:
    """
    Reverse translates a protein sequence into a DNA sequence and chooses an RBS using the highest CAI codon for each amino acid.
    """

    def __init__(self):
        self.codonUsage = {}
        self.rbsChooser = None

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
        # Translate peptide to codons
        codons = [self.aminoAcidToCodon[aa] for aa in peptide]

        # Append the stop codon (TAA in this case)
        codons.append("TAA")

        # Build the CDS from the codons
        cds = ''.join(codons)

        # Choose an RBS
        selectedRBS = self.rbsChooser.run(cds, ignores)

        # Return the Transcript object
        return Transcript(selectedRBS, peptide, codons)

if __name__ == "__main__":
    # Example usage of TranscriptDesigner
    peptide = "MYPFIRTARMTV"
    
    designer = TranscriptDesigner()
    designer.initiate()

    ignores = set()
    transcript = designer.run(peptide, ignores)
    
    # Print out the transcript information
    print(transcript)
