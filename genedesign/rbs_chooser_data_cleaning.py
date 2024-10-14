# Read and extract the relevant gene data
from Bio import SeqIO
from collections import defaultdict
import pandas as pd
# Function to extract UTR, gene, and CDS information from the GenBank file
def extract_genes_info(genbank_file):
    gene_dict = defaultdict(dict)  # Dictionary to store gene info
    for record in SeqIO.parse(genbank_file, "genbank"):
        for feature in record.features:
            if feature.type == "gene":
                locus_tag = feature.qualifiers.get("locus_tag", [None])[0]
                gene_name = feature.qualifiers.get("gene", [None])[0]

                # CDS information
                cds_feature = None
                for cds in record.features:
                    if cds.type == "CDS" and cds.qualifiers.get("locus_tag") == [locus_tag]:
                        cds_feature = cds
                        break

                if cds_feature:
                    start, end = cds_feature.location.start, cds_feature.location.end
                    strand = cds_feature.location.strand
                    if strand == 1:  # Forward strand
                        utr_start = max(0, start - 50)
                        utr_seq = record.seq[utr_start:start]
                    else:  # Reverse strand, we need to reverse complement
                        utr_start = end
                        utr_seq = record.seq[utr_start:utr_start + 50].reverse_complement()

                    cds_seq = cds_feature.extract(record.seq)
                    # Save the gene information in the dictionary
                    gene_dict[locus_tag] = {
                        "gene": gene_name,
                        "UTR": utr_seq,
                        "CDS": cds_seq
                    }
    return gene_dict

# Example usage
genbank_file = "genedesign/data/sequence.gb"  # Now using the correct file name
genes_info = extract_genes_info(genbank_file)