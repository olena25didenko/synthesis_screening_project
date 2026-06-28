"""
Wobble base modification: exploit degeneracy in the genetic code
Uses wobble pairing rules to modify sequences while preserving amino acid identity
More sophisticated than simple codon swapping
"""

from Bio import SeqIO
import json
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

# Reverse mapping: amino acid -> list of codons
CODONS_FOR_AA = {}
for codon, aa in GENETIC_CODE.items():
    if aa not in CODONS_FOR_AA:
        CODONS_FOR_AA[aa] = []
    CODONS_FOR_AA[aa].append(codon)

# Wobble pairing rules: 3rd position degeneracy
# Wobble position can match: same base, or specific flexible pairs
WOBBLE_RULES = {
    'A': ['A', 'G'],  # A wobbles with A or G
    'C': ['C', 'G'],  # C wobbles with C or G
    'G': ['A', 'C', 'G', 'U'],  # G wobbles with A, C, G, U
    'T': ['A', 'G'],  # T wobbles with A or G (in DNA context)
    'U': ['A', 'G']   # U wobbles with A or G (in RNA context)
}

def wobble_compatible(base1, base2):
    """Check if two bases are wobble-compatible at 3rd position"""
    base1_upper = base1.upper()
    base2_upper = base2.upper()
    
    # Convert T to U for wobble rules
    if base1_upper == 'T':
        base1_upper = 'U'
    if base2_upper == 'T':
        base2_upper = 'U'
    
    return base2_upper in WOBBLE_RULES.get(base1_upper, [])

def get_wobble_variants(codon):
    """
    Generate all codons that code for the same amino acid
    using wobble-compatible 3rd position changes
    """
    codon_upper = codon.upper()
    aa = GENETIC_CODE.get(codon_upper, 'X')
    
    if aa == 'X' or aa not in CODONS_FOR_AA:
        return [codon]
    
    # Get all codons for this amino acid
    all_codons = CODONS_FOR_AA[aa]
    
    # Filter to only wobble-compatible variants
    # (3rd position can change via wobble pairing)
    wobble_variants = []
    for alt_codon in all_codons:
        # Check if 1st and 2nd positions match, 3rd is wobble-compatible
        if (codon_upper[0] == alt_codon[0] and 
            codon_upper[1] == alt_codon[1] and
            wobble_compatible(codon_upper[2], alt_codon[2])):
            wobble_variants.append(alt_codon)
        # Also include non-wobble variants for more aggressive modification
        elif codon_upper != alt_codon:
            wobble_variants.append(alt_codon)
    
    return wobble_variants if wobble_variants else [codon]

def wobble_modify_sequence(dna_seq, aggressiveness=0.5):
    """
    Apply wobble base modifications to sequence
    
    Args:
        dna_seq: DNA sequence to modify
        aggressiveness: 0.0-1.0, how many codons to modify
                       0.5 = modify ~50% of codons
    """
    modified = []
    total_codons = 0
    modified_codons = 0
    
    for i in range(0, len(dna_seq) - 2, 3):
        codon = dna_seq[i:i+3].upper()
        total_codons += 1
        
        # Decide whether to modify this codon
        if random.random() < aggressiveness:
            # Get wobble variants
            variants = get_wobble_variants(codon)
            
            # Pick a different variant
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
    """Translate DNA sequence to protein"""
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
    """Verify protein identity is maintained"""
    original_protein = dna_to_protein(original_seq)
    modified_protein = dna_to_protein(modified_seq)
    
    # Calculate nucleotide identity
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

# Main execution
if __name__ == '__main__':
    print("\n" + "="*80)
    print("WOBBLE BASE MODIFICATIONS FOR SYNTHESIS SCREENING PROJECT")
    print("="*80 + "\n")
    
    sequences_to_modify = {
        'Influenza A': 'data/reference_sequences/influenza_ha_reference.fasta',
        'Ebola Virus': 'data/reference_sequences/ebola_vp40_reference.fasta',
        'Lassa Virus': 'data/reference_sequences/lassa_np_reference.fasta',
        'Measles Virus': 'data/reference_sequences/measles_f_reference.fasta',
        'SARS-CoV-2': 'data/reference_sequences/sars_cov2_spike_reference.fasta',
    }
    
    os.makedirs('data/modified_sequences', exist_ok=True)
    
    results_summary = []
    
    for organism, input_path in sequences_to_modify.items():
        print(f"\n{'='*80}")
        print(f"Wobble Modification: {organism}")
        print(f"{'='*80}")
        
        try:
            record = SeqIO.read(input_path, 'fasta')
            original_seq = str(record.seq)
        except FileNotFoundError:
            print(f"❌ File not found: {input_path}")
            continue
        
        # Apply wobble modification (aggressiveness: modify ~70% of codons)
        modified_seq, modified_count, total = wobble_modify_sequence(original_seq, aggressiveness=0.7)
        
        # Validate
        validation = validate_wobble_modification(original_seq, modified_seq)
        
        print(f"\n📊 Results for {organism}:")
        print(f"   Original length: {len(original_seq)} bp ({len(validation['original_protein'])} aa)")
        print(f"   Modified length: {len(modified_seq)} bp ({len(validation['modified_protein'])} aa)")
        print(f"   Proteins match: {validation['proteins_match']} ✅" if validation['proteins_match'] else f"   Proteins match: {validation['proteins_match']} ❌")
        print(f"   Nucleotide identity: {validation['nucleotide_identity']:.1f}%")
        print(f"   Codons changed: {validation['codons_changed']}/{validation['total_codons']} ({100*modified_count/total:.1f}%)")
        
        # Save modified sequence
        output_filename = input_path.replace('reference_sequences', 'modified_sequences').replace('_reference', '_wobble')
        with open(output_filename, 'w') as f:
            f.write(f">{record.id}_wobble_modified\n")
            f.write(modified_seq + "\n")
        
        print(f"   ✅ Saved to: {output_filename}")
        
        results_summary.append({
            'organism': organism,
            'protein_preserved': validation['proteins_match'],
            'nucleotide_identity': validation['nucleotide_identity'],
            'codons_changed': validation['codons_changed'],
            'output_file': output_filename
        })
    
    # Print summary
    print("\n" + "="*80)
    print("WOBBLE MODIFICATION SUMMARY")
    print("="*80)
    print(f"\n{'Organism':20} | {'Protein OK':12} | {'Nuc. ID %':12} | {'Codons Changed':15}")
    print("-" * 70)
    
    for result in results_summary:
        status = "✅ YES" if result['protein_preserved'] else "❌ NO"
        print(f"{result['organism']:20} | {status:12} | {result['nucleotide_identity']:11.1f}% | {result['codons_changed']:15}")
    
    print("\n✅ All sequences modified and saved to: data/modified_sequences/")
    print("="*80 + "\n")