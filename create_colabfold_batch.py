from Bio import SeqIO
import os
import json

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

# Create batch file for ColabFold with proper formatting
print('Creating ColabFold batch input file...')

# Original sequences
originals = {}
ref_dir = 'data/reference_sequences'
for fasta_file in sorted([f for f in os.listdir(ref_dir) if f.endswith('.fasta') and 'all_' not in f]):
    for record in SeqIO.parse(os.path.join(ref_dir, fasta_file), 'fasta'):
        protein = dna_to_protein(str(record.seq))
        base_name = fasta_file.replace('.fasta', '')
        originals[base_name] = (protein, record.id)

# Create FASTA with proper headers (ColabFold expects FASTA format with : for chains)
output = open('colabfold_batch_input.fasta', 'w')

print(f'\n{'='*80}')
print(f'STRUCTURAL VALIDATION BATCH INPUT')
print(f'{'='*80}\n')

count = 0

# Original sequences
print('Adding original sequences...')
for base_name in sorted(originals.keys()):
    protein, record_id = originals[base_name]
    output.write(f'>{base_name}_ORIGINAL\n{protein}\n')
    count += 1

# Modified sequences
print('Adding modified sequences...')
mod_dir = 'data/modified_sequences'
for fasta_file in sorted(os.listdir(mod_dir)):
    if fasta_file.endswith('.fasta'):
        for record in SeqIO.parse(os.path.join(mod_dir, fasta_file), 'fasta'):
            protein = dna_to_protein(str(record.seq))
            if protein:  # Only add if translation successful
                clean_id = fasta_file.replace('.fasta', '').upper()
                output.write(f'>{clean_id}\n{protein}\n')
                count += 1

output.close()

print(f'\n✅ Total sequences in batch: {count}')
print(f'✅ File: colabfold_batch_input.fasta')
print(f'\nNext steps:')
print(f'  1. Upload colabfold_batch_input.fasta to ColabFold')
print(f'  2. Use batch mode: Runtime → Run all')
print(f'  3. Download result ZIP with predictions + pLDDT scores')
print(f'  4. Extract and analyze structural conservation')

# Show file size
import os
size_kb = os.path.getsize('colabfold_batch_input.fasta') / 1024
print(f'\n✅ File size: {size_kb:.1f} KB')
