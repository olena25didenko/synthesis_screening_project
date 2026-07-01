from Bio import SeqIO
import os

GENETIC_CODE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}

def dna_to_protein(dna_seq):
    protein = []
    for i in range(0, len(dna_seq) - 2, 3):
        codon = dna_seq[i:i+3].upper()
        if codon in GENETIC_CODE:
            aa = GENETIC_CODE[codon]
            if aa == '*':
                break
            protein.append(aa)
        else:
            protein.append('X')
    return ''.join(protein)

# Get all sequences
all_sequences = []

# Originals
ref_dir = 'data/reference_sequences'
for fasta_file in sorted([f for f in os.listdir(ref_dir) if f.endswith('.fasta') and 'all_' not in f]):
    for record in SeqIO.parse(os.path.join(ref_dir, fasta_file), 'fasta'):
        protein = dna_to_protein(str(record.seq))
        all_sequences.append((f'{fasta_file.replace(\".fasta\", \"\")}_ORIGINAL', protein))

# Modified
mod_dir = 'data/modified_sequences'
for fasta_file in sorted(os.listdir(mod_dir)):
    if fasta_file.endswith('.fasta'):
        for record in SeqIO.parse(os.path.join(mod_dir, fasta_file), 'fasta'):
            protein = dna_to_protein(str(record.seq))
            all_sequences.append((fasta_file.replace('.fasta', '').upper(), protein))

# Create 4 batch files (50 sequences each)
batch_size = 50
num_batches = (len(all_sequences) + batch_size - 1) // batch_size

for batch_num in range(num_batches):
    start_idx = batch_num * batch_size
    end_idx = min(start_idx + batch_size, len(all_sequences))
    batch_sequences = all_sequences[start_idx:end_idx]
    
    filename = f'colabfold_batch_{batch_num + 1}.fasta'
    with open(filename, 'w') as f:
        for seq_id, protein in batch_sequences:
            f.write(f'>{seq_id}\n{protein}\n')
    
    print(f'✅ {filename}: {len(batch_sequences)} sequences')

print(f'\n✅ Created {num_batches} batch files')
print(f'✅ Total sequences: {len(all_sequences)}')
print(f'\nInstructions:')
print(f'  1. Create 4 new ColabFold tabs')
print(f'  2. For each batch file:')
print(f'     - Upload colabfold_batch_N.fasta')
print(f'     - Set jobname = \"batch_N\"')
print(f'     - Runtime → Run all')
print(f'     - Download result ZIP')
