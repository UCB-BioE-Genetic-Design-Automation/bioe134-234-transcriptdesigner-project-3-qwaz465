class InternalRBSChecker:
    def __init__(self):
        #add vars here
        self.sd_sequence = None
        self.start_codon = None
        self.sd_length = None
    
    def initiate(self):
        self.sd_sequence = 'AGGAGG'
        self.start_codon = 'ATG'
        self.sd_length = len(self.sd_sequence)

    def run(self, cds):
        # Search for SD sequence in the string
        sd_index = cds.find(self.sd_sequence)
        # If SD sequence is found, check for start codon 8 nucleotides later
        while sd_index != -1:
            # Calculate the position where start codon should be
            start_index = sd_index + self.sd_length + 8  # 6 nucleotides for SD, plus 8 spacer nucleotides

            # Check if there's a start codon at this position
            if cds[start_index:start_index + 3] == self.start_codon:
                return False, cds[start_index:start_index + 3]  # SD sequence with start codon found

            # Search for the next occurrence of the SD sequence
            sd_index = cds.find(self.sd_sequence, sd_index + 1)

        # Return True if no SD+start codon sequence is found
        return True, None