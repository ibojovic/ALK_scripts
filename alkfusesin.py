import random
import pysam

# ALK isoform 2 precursor sequence (RNA sequence)
alk_rna_sequence = (
    "CAGCCTTCCCTGGCTCCCTCCCCATTTCCTCTCATGGGCATTTCTTCTAATAAAATCTGCAGACCATA"
    "TTGGGTCTAATCCCATCTCCAGTCTGCTTCTTGGAGGAACCAGACTAACATGACTCTGCCCTATATAATA"
    "CAAATAATTATTTTCCATATATCTGATTTTTAGCTTTGCATTTACTTTAAATCATGCTTCAATTAATGACA"
    "CACCTTCTTTAATCATTTTATTAGTATTTCTAAGTATGATGGAAAGGTTCAGAGCTCAGGGGAGGATATGG"
    "AGATCCAGGGAGGCTTCCTGTAGGAAGTGGCCTGTGTAGTGCTTCAAGGGCCAGGCTGCCAGGCCATGTTG"
    "CAGCTGACCACCCACCTGCAGTGTACCGCCGGAAGCACCATGAGCTGCAAGCCATGCAGATGGAGCTGCAG"
    "AGCCCTGAGTACAAGCTGAGCAAGCTCCGCACCTCGACCATCATGACCGACTACAACCCAACTACTGCTTT"
    "GCTGGCAAGACCTCCTCCATCAGTGACCTGAAGGAGGTGCCGCGGAAAAACATCACCTCATTCCGGGTCTG"
    "GGCCATGGCGCCTTTGGGGAGGTGTATGAAGGCCAGGTGTCCGGAATGCCCAACGACCCAAGCCCCCTGCA"
    "AGTGGCTGTGAAGACGCTGCCTGAAGTGTGCTCTGAACAGGACGAACTGGATTTCCTCATGGAAGCCCTGA"
    "TCATCAGCAAATTCAACCACCAGAACATTGTTGCTGCATTGGGGTGAGCCTGCAATCCCTGCCCAGGTTCAT"
    "CCTGCTGGAGCTCATGGCGGGGGGAGACCTCAAGTCCTTCCTCCGAGAGACCCGCCCTCGCCCGAGCCAGCC"
    "CTCCTCCCTGGCCATGCTGGACCTTCTGCACGTGGCTCGGGACATTGCCTGTGGCTGTCAGTATTTGGAGGA"
    "AAACCACTTCATCCACCGAGACATTGCTGCCAGAAACTGCCTCTTGACCTGTCCAGGCCCTGGAAGAGTGGC"
    "CAAAGATTCGAGACTTCGGGATGGCCCGAGACATCTACAGGGCGAGCTACTATAGAAAGGGAGGCTGTGCCA"
    "TGCTGCCAGTTAAGTGGATGCCCCCAGAGGCCTTCATGGAAGGAATATTCACTTCTAAAAACAGACACATGG"
    "TCCTTTGGAGTGCTGCTATGGGAAATCTTTTCTCTTGGATATATGCCATACCCAGCAAAAGCAACCAGGAAG"
    "TTCTGGAGTTTGTCACCAGTGGAGGCCGGATGGACCCACCCAAGAACTGCCCTGGGCCTGTATACCGGATAA"
    "TGACTCAGTGCTGGCAACATCAGCCTGAAGACAGGCCCAACTTTGCCATCATTTTGGAGAGGATTGAATACT"
    "GCACCCAGGACCCGGATGTAATCAACACCGCTTTGCCGATAGAATATGGTCCACTTG"
)

# Function to create the fused RNA sequence of ALK with another gene
def create_fused_rna_sequence(alk_rna_sequence, other_gene_rna_sequence):
    # Choose a random fusion point in ALK and the other gene
    alk_fusion_point = random.randint(0, len(alk_rna_sequence))
    other_gene_fusion_point = random.randint(0, len(other_gene_rna_sequence))

    # Generate the fused RNA sequence by concatenating the portions before and after the fusion points
    fused_rna_sequence = alk_rna_sequence[:alk_fusion_point] + other_gene_rna_sequence[other_gene_fusion_point:]
    return fused_rna_sequence

# Function to retrieve gene RNA sequences based on gene names from a sequence database (e.g., NCBI GenBank)
def fetch_gene_sequences(gene_list_file):
    Entrez.email = "itanabojovic@gmail.com"  # Replace with your email address for Entrez API
    gene_rna_sequences = {}
    with open(gene_list_file, "r") as f:
        for line in f:
            gene_info = line.strip().split("\t")
            gene_name = gene_info[0]
            gene_id = gene_info[1] if len(gene_info) > 1 else None
            gene_rna_sequence = None
            try:
                # Fetch the gene RNA sequence from NCBI GenBank using gene name or gene ID
                handle = Entrez.efetch(db="nucleotide", id=gene_id, term=gene_name, rettype="gb", retmode="text")
                record = SeqIO.read(handle, "genbank")
                gene_rna_sequence = str(record.seq).replace("T", "U")  # Replace T with U for RNA
                handle.close()
            except Exception as e:
                print(f"Error fetching RNA sequence for {gene_name}: {str(e)}")
            if gene_rna_sequence:
                gene_rna_sequences[gene_name] = gene_rna_sequence
    return gene_rna_sequences

# Rest of the functions remain the same as before...

# Main function to generate ALK fusion RNA-seq FASTQ files for each gene in the list
def generate_fusion_fastq_files(gene_list_file, read_length, num_reads_per_gene):
    gene_rna_sequences = fetch_gene_sequences(gene_list_file)
    for gene_name, gene_rna_sequence in gene_rna_sequences.items():
        # Create the fused RNA sequence for the current gene
        fused_rna_sequence = create_fused_rna_sequence(alk_rna_sequence, gene_rna_sequence)

        # Simulate paired-end reads from the fused RNA sequence
        paired_reads = simulate_paired_end_reads(fused_rna_sequence, read_length, num_reads_per_gene)

        # Generate quality scores for simulated reads
        quality_scores = [(generate_quality_scores(read1), generate_quality_scores(read2)) for read1, read2 in paired_reads]

        # Output FASTQ file paths
        output_fastq_file1 = f"{gene_name}_ALK_fusion_R1.fastq"
        output_fastq_file2 = f"{gene_name}_ALK_fusion_R2.fastq"

        # Save the simulated paired-end reads and quality scores to FASTQ files
        save_to_fastq(output_fastq_file1, [read[0] for read in paired_reads], [quality[0] for quality in quality_scores])
        save_to_fastq(output_fastq_file2, [read[1] for read in paired_reads], [quality[1] for quality in quality_scores])

# Example usage with a gene list file with gene names (and optional gene IDs)
gene_list_file = "gene_list.txt"
# Format of gene_list.txt:
# GeneX	GeneX_ID
# GeneY	GeneY_ID

# Parameters for the simulation
read_length = 75
num_reads_per_gene = 1000

# Generate ALK fusion RNA-seq FASTQ files for each gene in the list
generate_fusion_fastq_files(gene_list_file, read_length, num_reads_per_gene)
