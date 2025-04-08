from flask import Flask, render_template, request
from needleman_wunsch import global_alignment
from smith_waterman import local_alignment
from gene_labeling import label_gene_regions
from graph_module import Graph, bfs_locus_to_disease
import json

app = Flask(__name__)

# Load reference wild-type sequence
with open("original_mtDNA/original_mtdna.fasta", "r") as f:
    wild_seq = "".join(line.strip() for line in f if not line.startswith(">"))

print(f"✅ Wild-type sequence loaded. First 50 bases: {wild_seq[:50]}...\n")

# Load reference gene locus positions
with open("mt_locus.json") as f:
    gene_locus = json.load(f)

# Initialize the Graph for disease prediction
graph = Graph()
graph.load_from_csv("mito.csv")  # Ensure mito.csv is present in the same directory

@app.route("/", methods=["GET", "POST"])
def home():
    if request.method == "POST":
        file = request.files.get("usersequence")
        if not file:
            return "Error: No file uploaded", 400

        mutated_seq = file.read().decode("utf-8").strip()
        mutated_seq = "".join(line.strip() for line in mutated_seq.split("\n") if not line.startswith(">"))

        print(f"✅ Mutated sequence uploaded. First 50 bases: {mutated_seq[:50]}...\n")
        print("🔄 Performing Needleman-Wunsch Global Alignment...\n")

        # Perform global alignment
        wild_aligned, mutated_aligned, global_score = global_alignment(wild_seq, mutated_seq)
        print(f"✅ Global Alignment Complete. Score: {global_score}\n")

        # Identify gene regions
        labeled_sequences = label_gene_regions(wild_aligned, mutated_aligned)

        # Perform local alignment (Smith-Waterman) for all loci
        local_results = {}
        mutation_impact = {} 
        all_mutations = []

        print("🔄 Performing Smith-Waterman Local Alignment for all loci...\n")
        for locus, sequences in labeled_sequences.items():
            wild_locus_seq = sequences["wild_seq"]
            mutated_locus_seq = sequences["mutated_seq"]

            local_aligned_seq1, local_aligned_seq2, local_score = local_alignment(wild_locus_seq, mutated_locus_seq)

            # Detect mutations and format them correctly
            mutations = []
            for i, (a, b) in enumerate(zip(local_aligned_seq1, local_aligned_seq2)):
                if a != b:
                    position = gene_locus[locus][0] + i  # Adjust position to match genome
                    if a == "-":
                        mutation = f"m.{position}_{position+1}ins{b}"
                    elif b == "-":
                        mutation = f"m.{position}{a}_del"
                    else:
                        mutation = f"m.{position}{a}>{b}"
                    mutations.append(mutation)
                    all_mutations.append(mutation)

            num_mutations = len(mutations)
            mutation_impact[locus] = {"score": local_score, "mutations": num_mutations, "alleles": mutations}

            local_results[locus] = {
                "score": local_score,
                "aligned_seq1": local_aligned_seq1,
                "aligned_seq2": local_aligned_seq2,
                "mutations": mutations
            }

            print(f"✅ Local Alignment Complete for {locus}. Score: {local_score}, Mutations: {num_mutations}")
            print(f"🧬 Detected Mutations: {mutations}\n")

        # Predict diseases for ALL mutations (no top allele nonsense)
        predicted_diseases = set()

        print("🔍 Checking all mutations in graph...\n")
        for mutation in all_mutations:
            if mutation in graph.graph:
                print(f"✅ {mutation} found in graph. Running BFS...")
                paths = bfs_locus_to_disease(graph.graph, mutation)
                for path in paths:
                    if len(path) > 1:
                        predicted_diseases.add(path[-1])  # disease at the end
            else:
                print(f"❌ {mutation} not found in graph.")

        print(f"⚠️ Final Predicted Diseases: {predicted_diseases}")

        return render_template("result.html",
                               global_score=global_score,
                               aligned_seq1=wild_aligned,
                               aligned_seq2=mutated_aligned,
                               local_results=local_results,
                               mutation_impact=mutation_impact,
                               predicted_diseases=list(predicted_diseases),
                               all_mutations=all_mutations)

    return render_template("index.html")

if __name__ == "__main__":
    print("🚀 Flask server is running at http://127.0.0.1:5000/\n")
    app.run(debug=True)