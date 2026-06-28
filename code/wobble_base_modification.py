"""
Wobble base modifications - AUTO-DISCOVERS all FASTA files
"""

from Bio import SeqIO
import os
import random

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

CODONS_FOR_AA = {}
for codon, aa in GENETIC_CODE.items():
    if aa not in CODONS_FOR_AA:
        CODONS_FOR_AA[aa] = []
    CODONS_FOR_AA[aa].append(codon)

WOBBLE_RULES = {
    'A': ['A', 'G'], 'C': ['C', 'G'], 'G': ['A', 'C', 'G', 'U'],
    'T': ['A', 'G'], 'U': ['A', 'G']
}

def wobble_compatible(base1, base2):
    base1_upper = base1.upper()
    base2_upper = base2.upper()
    if base1_upper == 'T':
        base1_upper = 'U'
    if base2_upper == 'T':
        base2_upper = 'U'
    return base2_upper in WOBBLE_RULES.get(base1_upper, [])

def get_wobble_variants(codon):
    codon_upper = codon.upper()
    aa = GENETIC_CODE.get(codon_upper, 'X')
    if aa == 'X' or aa not in CODONS_FOR_AA:
        return [codon]
    all_codons = CODONS_FOR_AA[aa]
    wobble_variants = []
    for alt_codon in all_codons:
        if (codon_upper[0] == alt_codon[0] and 
            codon_upper[1] == alt_codon[1] and
            wobble_compatible(codon_upper[2], alt_codon[2])):
            wobble_variants.append(alt_codon)
        elif codon_upper != alt_codon:
            wobble_variants.append(alt_codon)
    return wobble_variants if wobble_variants else [codon]

def wobble_modify_sequence(dna_seq, aggressiveness=0.7):
    modified = []
    total_codons = 0
    modified_codons = 0
    for i in range(0, len(dna_seq) - 2, 3):
        codon = dna_seq[i:i+3].upper()
        total_codons += 1
        if random.random() < aggressiveness:
            variants = get_wobble_variants(codon)
            alternative_codons = [c for c in variants if c != codon]
            if alternative_codons:
                new_codon = random.choice(alternative_codons)
                modified.append(new_codon)
                modified_codons += 1
            else:
                modified.append(codon)
        else:
            modified.append(codon)
    return ''.join(modified), modified_codons, total_codons

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

def validate_wobble_modification(original_seq, modified_seq):
    original_protein = dna_to_protein(original_seq)
    modified_protein = dna_to_protein(modified_seq)
    matching_codons = sum(
        1 for i in range(0, min(len(original_seq), len(modified_seq)), 3)
        if original_seq[i:i+3].upper() == modified_seq[i:i+3].upper()
    )
    total_codons = min(len(original_seq), len(modified_seq)) // 3
    nucleotide_identity = (matching_codons / total_codons * 100) if total_codons > 0 else 0
    return {
        'proteins_match': original_protein == modified_protein,
        'original_protein': original_protein,
        'modified_protein': modified_protein,
        'nucleotide_identity': nucleotide_identity,
        'codons_changed': total_codons - matching_codons,
        'total_codons': total_codons
    }

if __name__ == '__main__':
    print("\n" + "="*80)
    print("WOBBLE BASE MODIFICATIONS")
    print("AUTO-DISCOVERING ALL SEQUENCES")
    print("="*80 + "\n")
    
    os.makedirs('data/modified_sequences', exist_ok=True)
    
    # AUTO-DISCOVER all FASTA files
    ref_dir = 'data/reference_sequences'
    fasta_files = [f for f in os.listdir(ref_dir) if f.endswith('.fasta') and 'all_sequences' not in f]
    
    print(f"Found {len(fasta_files)} FASTA files to process:\n")
    
    results_summary = []
    
    for fasta_file in sorted(fasta_files):
        input_path = os.path.join(ref_dir, fasta_file)
        
        print(f"Processing: {fasta_file}...", end=' ')
        
        try:
            record = SeqIO.read(input_path, 'fasta')
            original_seq = str(record.seq)
        except FileNotFoundError:
            print(f"❌ File not found")
            continue
        
        modified_seq, modified_count, total = wobble_modify_sequence(original_seq, aggressiveness=0.7)
        validation = validate_wobble_modification(original_seq, modified_seq)
        
        output_filename = fasta_file.replace('.fasta', '_wobble.fasta')
        output_path = os.path.join('data/modified_sequences', output_filename)
        
        with open(output_path, 'w') as f:
            f.write(f">{record.id}_wobble_modified\n")
            f.write(modified_seq + "\n")
        
        print(f"✅ ({validation['nucleotide_identity']:.1f}% identity)")
        
        results_summary.append({
            'file': fasta_file,
            'protein_preserved': validation['proteins_match'],
            'nucleotide_identity': validation['nucleotide_identity'],
            'codons_changed': validation['codons_changed'],
            'output_file': output_path
        })
    
    print("\n" + "="*80)
    print(f"WOBBLE MODIFICATION SUMMARY ({len(results_summary)} sequences)")
    print("="*80)
    print(f"\n{'File':40} | {'Protein OK':12} | {'Nuc. ID %':12}")
    print("-" * 70)
    
    for result in results_summary:
        status = "✅ YES" if result['protein_preserved'] else "❌ NO"
        print(f"{result['file']:40} | {status:12} | {result['nucleotide_identity']:11.1f}%")
    
    print(f"\n✅ Processed {len(results_summary)} sequences")
    print("="*80 + "\n")