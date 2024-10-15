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

def proteomics_prune(file_path):
    proteomics_data = pd.read_csv(file_path, sep="\t", names=["string_external_id", "abundance"])

    # Convert the abundance column to numeric, coerce errors (for any non-numeric values)
    proteomics_data['abundance'] = pd.to_numeric(proteomics_data['abundance'], errors='coerce')

    # Drop any rows with missing or NaN values in the abundance column
    proteomics_data.dropna(subset=['abundance'], inplace=True)

    # Extract the locus tag (string after the period) from string_external_id
    proteomics_data['locus_tag'] = proteomics_data['string_external_id'].apply(lambda x: x.split('.')[1])

    # Sort data by abundance in descending order
    proteomics_data = proteomics_data.sort_values(by='abundance', ascending=False)

    # Calculate the number of rows corresponding to the top 5%
    top_5_percent_count = int(0.05 * len(proteomics_data))

    # Select the top 5% most abundant locus tags
    top_5_percent = proteomics_data.head(top_5_percent_count)

    # Return as a list of tuples (locus tag, abundance)
    top_5_percent_list = list(top_5_percent[['locus_tag', 'abundance']].itertuples(index=False, name=None))
    return top_5_percent_list


# Example usage to be used in rbs_chooser ?
# genbank_file = "genedesign/data/sequence.gb"  # Now using the correct file name
# genes_info = extract_genes_info(genbank_file)
# file_path = 'genedesign/data/511145-WHOLE_ORGANISM-integrated.txt'
# top_5_percent_list = proteomics_prune(file_path)