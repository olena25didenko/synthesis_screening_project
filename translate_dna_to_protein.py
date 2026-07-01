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

# Create TWO files:
# 1. FASTA format (for documentation)
fasta_out = open('colabfold_proteins_batch.fasta', 'w')
# 2. Sequences only (for ColabFold)
seqs_out = open('colabfold_sequences_only.txt', 'w')

print('Translating DNA → Protein sequences...\n')

# Original sequences
ref_dir = 'data/reference_sequences'
count = 0
for fasta_file in sorted([f for f in os.listdir(ref_dir) if f.endswith('.fasta') and 'all_' not in f]):
    for record in SeqIO.parse(os.path.join(ref_dir, fasta_file), 'fasta'):
        protein = dna_to_protein(str(record.seq))
        fasta_out.write(f'>{record.id}\n{protein}\n')
        seqs_out.write(f'{protein}\n')
        count += 1

# Modified sequences
mod_dir = 'data/modified_sequences'
for fasta_file in sorted(os.listdir(mod_dir)):
    if fasta_file.endswith('.fasta'):
        for record in SeqIO.parse(os.path.join(mod_dir, fasta_file), 'fasta'):
            protein = dna_to_protein(str(record.seq))
            fasta_out.write(f'>{record.id}\n{protein}\n')
            seqs_out.write(f'{protein}\n')
            count += 1

fasta_out.close()
seqs_out.close()

print(f'✅ Translated {count} sequences')
print(f'✅ colabfold_proteins_batch.fasta (FASTA format)')
print(f'✅ colabfold_sequences_only.txt (sequences only for ColabFold)')
