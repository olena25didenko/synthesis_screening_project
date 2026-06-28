import os
from Bio import SeqIO

output = open('data/reference_sequences/all_49_sequences.fasta', 'w')
ref_dir = 'data/reference_sequences'
fasta_files = sorted([f for f in os.listdir(ref_dir) if f.endswith('.fasta') and 'all_' not in f])

count = 0
for fasta_file in fasta_files:
    for record in SeqIO.parse(os.path.join(ref_dir, fasta_file), 'fasta'):
        output.write(f'>{record.id}\n{record.seq}\n')
        count += 1

output.close()
print(f'✅ Combined {count} sequences from {len(fasta_files)} files')
