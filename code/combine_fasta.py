import os

# Directory with individual sequences
ref_dir = "data/reference_sequences/"
output_file = "data/reference_sequences/all_sequences.fasta"

# Combine all FASTA files
with open(output_file, 'w') as outf:
    for filename in sorted(os.listdir(ref_dir)):
        if filename.endswith(".fasta") and filename != "all_sequences.fasta":
            filepath = os.path.join(ref_dir, filename)
            with open(filepath, 'r') as inf:
                outf.write(inf.read())
                outf.write("\n")  # Add newline between files

print(f"Combined {len([f for f in os.listdir(ref_dir) if f.endswith('.fasta') and f != 'all_sequences.fasta'])} FASTA files into {output_file}")