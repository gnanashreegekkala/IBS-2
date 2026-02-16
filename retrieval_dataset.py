from Bio import Entrez, SeqIO

Entrez.email = "wby1606@gmail.com"

viruses = [
    "Variola virus",
    "Monkeypox virus",
    "Cowpox virus",
    "Vaccinia virus"
]

sequences_per_virus = 250

for virus in viruses:
    print(f"\nSearching for {virus}...")

    search_handle = Entrez.esearch(
        db="nucleotide",
        term=f"{virus}[Organism]",
        retmax=sequences_per_virus
    )
    search_results = Entrez.read(search_handle)
    ids = search_results["IdList"]

    print(f"Found {len(ids)} sequences")

    if len(ids) == 0:
        print("No sequences found, skipping.")
        continue

    fetch_handle = Entrez.efetch(
        db="nucleotide",
        id=ids,
        rettype="fasta",
        retmode="text"
    )

    sequences = list(SeqIO.parse(fetch_handle, "fasta"))

    filename = virus.replace(" ", "_") + ".fasta"
    SeqIO.write(sequences, filename, "fasta")

    print(f"Saved to {filename}")
