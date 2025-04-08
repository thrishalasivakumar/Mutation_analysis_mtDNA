import json

# Load mitochondrial gene reference table from JSON
with open("mt_locus.json") as f:
    gene_locus = json.load(f)

def label_gene_regions(wild_aligned, mutated_aligned):
    labeled_sequences = {}

    # Build position mapping from reference (ungapped) to aligned wild-type
    genome_pos = 0
    alignment_index_map = {}  # Maps actual mtDNA positions to aligned sequence index

    for i, base in enumerate(wild_aligned):
        if base != "-":
            genome_pos += 1
            alignment_index_map[genome_pos] = i  # mtDNA position => index in aligned sequence

    for gene, (start, end) in gene_locus.items():
        # Get the correct aligned slice for this gene
        if start in alignment_index_map and end in alignment_index_map:
            aligned_start = alignment_index_map[start]
            aligned_end = alignment_index_map[end] + 1  # +1 to include the end position

            wild_gene_seq = wild_aligned[aligned_start:aligned_end]
            mutated_gene_seq = mutated_aligned[aligned_start:aligned_end]

            labeled_sequences[gene] = {
                "wild_seq": wild_gene_seq,
                "mutated_seq": mutated_gene_seq,
                "gaps": mutated_gene_seq.count("-")
            }
        else:
            print(f"⚠️ Skipping gene {gene} due to missing alignment positions.")

    return labeled_sequences