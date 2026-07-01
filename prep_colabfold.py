import os
from Bio import SeqIO

# Output file
output = open('colabfold_validation_batch.fasta', 'w')

# Original sequences
ref_dir = 'data/reference_sequences'
for fasta_file in sorted([f for f in os.listdir(ref_dir) if f.endswith('.fasta') and 'all_' not in f]):
    for record in SeqIO.parse(os.path.join(ref_dir, fasta_file), 'fasta'):
        output.write(f'>{record.id}\n{record.seq}\n')

# Modified sequences
mod_dir = 'data/modified_sequences'
for fasta_file in sorted(os.listdir(mod_dir)):
    if fasta_file.endswith('.fasta'):
        for record in SeqIO.parse(os.path.join(mod_dir, fasta_file), 'fasta'):
            output.write(f'>{record.id}\n{record.seq}\n')

output.close()
print(f'✅ Created colabfold_validation_batch.fasta')
