from Bio import SeqIO

input_files = {
    "variola": "Variola_final.fasta",
    "monkeypox": "Monkeypox_final.fasta",
    "cowpox": "Cowpox_final.fasta",
    "vaccinia": "Vaccinia_final.fasta"
}

num_sequences = 5
fragment_length = 500

for virus, file in input_files.items():
    output_file = virus + "_blast_subset.fasta"
    records = list(SeqIO.parse(file, "fasta"))
    
    with open(output_file, "w") as out:
        for record in records[:num_sequences]:
            fragment_seq = str(record.seq)[:fragment_length]
            out.write(f">{record.id}\n{fragment_seq}\n")
    
    print(f"{virus}: {num_sequences} sequences of {fragment_length} bases saved to {output_file}")
