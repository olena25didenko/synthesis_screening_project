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

# Spot-check: Just the 5 wobble-detected sequences + originals
validation_set = [
    'data/reference_sequences/hcov_nl63_spike.fasta',
    'data/reference_sequences/sars_cov1_spike_ay291451.fasta',
    'data/reference_sequences/mers_spike_kc881005.fasta',
    'data/reference_sequences/sars_cov2_spike_ov109917.fasta',
    'data/reference_sequences/marburg_gp.fasta',
    'data/modified_sequences/hcov_nl63_spike_wobble.fasta',
    'data/modified_sequences/sars_cov1_spike_ay291451_wobble.fasta',
    'data/modified_sequences/mers_spike_kc881005_wobble.fasta',
    'data/modified_sequences/sars_cov2_spike_ov109917_wobble.fasta',
    'data/modified_sequences/marburg_gp_wobble.fasta',
]

# Create validation FASTA
output = open('colabfold_spot_check.fasta', 'w')

for fasta_path in validation_set:
    if os.path.exists(fasta_path):
        for record in SeqIO.parse(fasta_path, 'fasta'):
            protein = dna_to_protein(str(record.seq))
            output.write(f'>{record.id}\n{protein}\n')

output.close()
print('✅ Created colabfold_spot_check.fasta')
print('✅ Contains: 5 originals + 5 wobble modifications')
print('✅ Ready for ColabFold upload')
