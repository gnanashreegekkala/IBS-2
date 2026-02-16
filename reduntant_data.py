from Bio import SeqIO
import numpy as np
from sklearn.metrics.pairwise import cosine_similarity
from itertools import product
import os

input_files = {
    "variola": "variola_virus.fasta",
    "monkeypox": "monkeypox_virus.fasta",
    "cowpox": "cowpox_virus.fasta",
    "vaccinia": "vaccinia_virus.fasta"
}

k = 3
similarity_threshold = 0.98

def generate_kmers(k):
    return [''.join(p) for p in product('ATGC', repeat=k)]

all_kmers = generate_kmers(k)
kmer_index = {kmer: idx for idx, kmer in enumerate(all_kmers)}

def seq_to_vector(seq):
    vector = np.zeros(len(all_kmers))
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i+k]
        if kmer in kmer_index:
            vector[kmer_index[kmer]] += 1
    return vector

for virus_name, file_path in input_files.items():
    
    records = list(SeqIO.parse(file_path, "fasta"))
    vectors = []
    kept_records = []

    for record in records:
        vec = seq_to_vector(str(record.seq).upper())
        
        if len(vectors) == 0:
            vectors.append(vec)
            kept_records.append(record)
        else:
            sims = cosine_similarity([vec], vectors)[0]
            if max(sims) < similarity_threshold:
                vectors.append(vec)
                kept_records.append(record)

    print("------", virus_name.upper(), "------")
    print("Original:", len(records))
    print("Non-redundant:", len(kept_records))
    print()

    output_file = virus_name + "_nonredundant.fasta"
    SeqIO.write(kept_records, output_file, "fasta")
