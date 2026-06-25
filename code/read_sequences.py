from Bio import SeqIO
import os

# Directory with reference sequences
ref_dir = "data/reference_sequences/"

print("=== Reading Reference Sequences ===\n")

for filename in os.listdir(ref_dir):
    if filename.endswith(".fasta"):
        filepath = os.path.join(ref_dir, filename)
        
        # Read FASTA file
        record = SeqIO.read(filepath, "fasta")
        
        # Print basic info
        print(f"File: {filename}")
        print(f"  Sequence ID: {record.id}")
        print(f"  Description: {record.description}")
        print(f"  Length: {len(record.seq)} bp")
        print(f"  First 50 bp: {str(record.seq[:50])}")
        print()