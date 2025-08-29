# Step 1: Download representative AChE protein sequences from NCBI into a FASTA file

from Bio import Entrez, SeqIO
import time

# 0) REQUIRED: tell NCBI who you are (replace with your details)
Entrez.email = "muskan0072003@gmail.com"     # <-- change this
Entrez.api_key = None                       # optional: "YOUR_NCBI_API_KEY"

# 1) Species and search terms (gene names differ a bit across species)
species_terms = {
    "Homo sapiens": "(ACHE[Gene] OR acetylcholinesterase[Protein Name])",
    "Mus musculus": "(Ache[Gene] OR acetylcholinesterase[Protein Name])",
    "Rattus norvegicus": "(Ache[Gene] OR acetylcholinesterase[Protein Name])",
    "Danio rerio": "(ache[Gene] OR acetylcholinesterase[Protein Name])",
    "Drosophila melanogaster": "(Ace[Gene] OR acetylcholinesterase[Protein Name])",
}

def search_and_fetch_sequences(organism, query_core, retmax=50):
    """Search NCBI Protein for AChE in a given organism, then fetch FASTA for the hits."""
    term = f'({query_core}) AND "{organism}"[Organism] AND refseq[filter]'
    with Entrez.esearch(db="protein", term=term, retmax=retmax) as h:
        result = Entrez.read(h)
    ids = result.get("IdList", [])
    if not ids:
        return []

    # Fetch all in one go (more polite to NCBI + faster)
    time.sleep(0.34)  # be nice to NCBI; with API key you can reduce this
    with Entrez.efetch(db="protein", id=",".join(ids), rettype="fasta", retmode="text") as h:
        records = list(SeqIO.parse(h, "fasta"))
    return records

def pick_representative(records):
    """
    Choose one good sequence:
    - Prefer curated RefSeq (NP_*)
    - Avoid 'partial' sequences
    - If multiple remain, pick the longest (usually canonical)
    """
    curated = [r for r in records if r.id.startswith("NP_")]  # curated RefSeq
    candidates = curated if curated else records
    candidates = [r for r in candidates if "partial" not in r.description.lower()] or candidates
    return max(candidates, key=lambda r: len(r.seq))

out = []
for org, core in species_terms.items():
    recs = search_and_fetch_sequences(org, core)
    if not recs:
        print(f"[WARN] No sequences found for {org}")
        continue
    best = pick_representative(recs)
    acc = best.id
    # Make a clean header like: Homo_sapiens|NP_000000.1|len=614
    best.id = f"{org.replace(' ', '_')}|{acc}|len={len(best.seq)}"
    best.description = ""
    out.append(best)
    print(f"Selected {org}: {acc} (length {len(best.seq)})")
    time.sleep(0.34)

# Save to a single multi-FASTA file
if out:
    SeqIO.write(out, "ache_proteins.fasta", "fasta")
    print(f"Saved {len(out)} sequences to ache_proteins.fasta")
else:
    print("No sequences selected.") 



#Step 2 — Align the AChE protein sequences

from Bio.Align.Applications import ClustalOmegaCommandline
from Bio import AlignIO

# Input FASTA (from Step 1)
in_file = "ache_proteins.fasta"
# Output files
aln_file = "ache_proteins.aln"
phy_file = "ache_proteins.phy"

# Run Clustal Omega
clustalomega_cline = ClustalOmegaCommandline(
    infile=in_file,
    outfile=aln_file,
    verbose=True,
    auto=True,
    force=True,        # overwrite if file exists
    outfmt="clu"       # alignment format: clustal
)
stdout, stderr = clustalomega_cline()

print("Alignment complete! Saved to:", aln_file)

# Convert to PHYLIP format (for phylogenetic tree later)
alignment = AlignIO.read(aln_file, "clustal")
AlignIO.write(alignment, phy_file, "phylip")
print("Also saved PHYLIP format:", phy_file) 



# Step 3 Find conserved regions in AChE Alignment

from Bio import AlignIO

# Load the alignment file
alignment = AlignIO.read("ache_proteins.aln", "clustal")

num_sequences = len(alignment)
alignment_length = alignment.get_alignment_length()

print(f"Number of sequences: {num_sequences}")
print(f"Alignment length: {alignment_length}")

# Step a: Find conserved columns
conserved_positions = []
for i in range(alignment_length):
    column = alignment[:, i]  # all residues at position i
    if "-" not in column:  # skip gaps
        if len(set(column)) == 1:  # all residues same 
            conserved_positions.append(i)

print(f"Total conserved positions: {len(conserved_positions)}")

# Step b: Extract conserved motifs (continuous blocks of conserved residues)
motifs = []
current = []
for pos in conserved_positions:
    if current and pos == current[-1] + 1:
        current.append(pos)
    else:
        if current:
            motifs.append(current)
        current = [pos]
if current:
    motifs.append(current)

print("\nConserved motifs found:")
for m in motifs:
    start, end = m[0]+1, m[-1]+1  # +1 for human-readable indexing 
    seq = alignment[0, m[0]:m[-1]+1]  # use first sequence as reference
    print(f"Positions {start}-{end}: {seq}")  


# STEP 4: Build and save a phylogenetic tree from the AChE alignment

from Bio import AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import matplotlib.pyplot as plt

# 1) Load your alignment (from Step 2)
aln_path = "ache_proteins.aln"     # clustal format
alignment = AlignIO.read(aln_path, "clustal")

print(f"Loaded alignment with {len(alignment)} sequences and length {alignment.get_alignment_length()} aa")

# 2) Compute pairwise distances
# For proteins, 'blosum62' is a good default; you can also try 'identity' to sanity-check.
calculator = DistanceCalculator('blosum62')
dm = calculator.get_distance(alignment)
print("\nPairwise distance matrix:")
print(dm)

# 3) Build trees
constructor = DistanceTreeConstructor()
nj_tree = constructor.nj(dm)        # Neighbor-Joining
upgma_tree = constructor.upgma(dm)  # Optional: UPGMA for comparison

# 4) Save trees to Newick (good to keep in your repo)
Phylo.write(nj_tree, "ache_tree_NJ.newick", "newick")
Phylo.write(upgma_tree, "ache_tree_UPGMA.newick", "newick")
print("\nSaved trees: ache_tree_NJ.newick, ache_tree_UPGMA.newick")

# 5) Render and save a figure (PNG)
plt.figure(figsize=(8, 5))
Phylo.draw(nj_tree, do_show=False)         # draws on current matplotlib axes
plt.tight_layout()
plt.savefig("ache_tree_NJ.png", dpi=300)
plt.close()
print("Saved figure: ache_tree_NJ.png")

# (Optional) also save a PDF
plt.figure(figsize=(8, 5))
Phylo.draw(nj_tree, do_show=False)
plt.tight_layout()
plt.savefig("ache_tree_NJ.pdf")
plt.close()
print("Saved figure: ache_tree_NJ.pdf")

print("\nDone! Use the NJ tree for your report; UPGMA is just for comparison.")


