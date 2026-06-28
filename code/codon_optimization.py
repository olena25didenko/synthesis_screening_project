"""
Codon optimization for synthesis screening
Uses the optimization maps to replace common codons with rare alternatives
"""

from Bio import SeqIO
import json
import os

# Standard genetic code
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

def optimize_sequence(dna_seq, optimization_map):
    """
    Optimize codon usage using the optimization map
    Replace common codons with rare alternatives
    """
    optimized = []
    
    for i in range(0, len(dna_seq) - 2, 3):
        codon = dna_seq[i:i+3].upper()
        
        # If this codon is in the optimization map, replace it
        if codon in optimization_map:
            optimized.append(optimization_map[codon])
        else:
            # Keep original if not in map
            optimized.append(codon)
    
    return ''.join(optimized)

def validate_optimization(original_seq, optimized_seq):
    """Verify protein identity is maintained"""
    original_protein = dna_to_protein(original_seq)
    optimized_protein = dna_to_protein(optimized_seq)
    
    # Calculate nucleotide identity
    matching_codons = sum(
        1 for i in range(0, min(len(original_seq), len(optimized_seq)), 3)
        if original_seq[i:i+3].upper() == optimized_seq[i:i+3].upper()
    )
    total_codons = min(len(original_seq), len(optimized_seq)) // 3
    
    nucleotide_identity = (matching_codons / total_codons * 100) if total_codons > 0 else 0
    
    return {
        'proteins_match': original_protein == optimized_protein,
        'original_protein': original_protein,
        'optimized_protein': optimized_protein,
        'nucleotide_identity': nucleotide_identity,
        'codons_changed': total_codons - matching_codons,
        'total_codons': total_codons
    }

# Main execution
if __name__ == '__main__':
    print("\n" + "="*80)
    print("CODON OPTIMIZATION FOR SYNTHESIS SCREENING PROJECT")
    print("="*80 + "\n")
    
    # Load optimization maps
    try:
        with open('data/codon_usage/optimization_maps.json', 'r') as f:
            all_maps = json.load(f)
    except FileNotFoundError:
        print("❌ optimization_maps.json not found. Run extract_codon_usage.py first!")
        exit(1)
    
    # Define which sequences to optimize
    sequences_to_optimize = {
        'Influenza A': ('data/reference_sequences/influenza_ha_reference.fasta', 'influenza_ha_optimized.fasta'),
        'Ebola Virus': ('data/reference_sequences/ebola_vp40_reference.fasta', 'ebola_vp40_optimized.fasta'),
        'Lassa Virus': ('data/reference_sequences/lassa_np_reference.fasta', 'lassa_np_optimized.fasta'),
        'Measles Virus': ('data/reference_sequences/measles_f_reference.fasta', 'measles_f_optimized.fasta'),
        'SARS-CoV-2': ('data/reference_sequences/sars_cov2_spike_reference.fasta', 'sars_cov2_spike_optimized.fasta'),
    }
    
    # Create results directory if it doesn't exist
    os.makedirs('data/modified_sequences', exist_ok=True)
    
    # Process each sequence
    results_summary = []
    
    for organism, (input_path, output_filename) in sequences_to_optimize.items():
        print(f"\n{'='*80}")
        print(f"Optimizing: {organism}")
        print(f"{'='*80}")
        
        # Check if optimization map exists
        if organism not in all_maps:
            print(f"⚠️  No optimization map for {organism}, skipping...")
            continue
        
        # Read original sequence
        try:
            record = SeqIO.read(input_path, 'fasta')
            original_seq = str(record.seq)
        except FileNotFoundError:
            print(f"❌ Sequence file not found: {input_path}")
            continue
        
        # Optimize
        optimization_map = all_maps[organism]
        optimized_seq = optimize_sequence(original_seq, optimization_map)
        
        # Validate
        validation = validate_optimization(original_seq, optimized_seq)
        
        print(f"\n📊 Results for {organism}:")
        print(f"   Original length: {len(original_seq)} bp ({len(validation['original_protein'])} aa)")
        print(f"   Optimized length: {len(optimized_seq)} bp ({len(validation['optimized_protein'])} aa)")
        print(f"   Proteins match: {validation['proteins_match']} ✅" if validation['proteins_match'] else f"   Proteins match: {validation['proteins_match']} ❌")
        print(f"   Nucleotide identity: {validation['nucleotide_identity']:.1f}%")
        print(f"   Codons changed: {validation['codons_changed']}/{validation['total_codons']}")
        
        # Save optimized sequence
        output_path = os.path.join('data/modified_sequences', output_filename)
        with open(output_path, 'w') as f:
            f.write(f">{record.id}_codon_optimized\n")
            f.write(optimized_seq + "\n")
        
        print(f"   ✅ Saved to: {output_path}")
        
        # Store result
        results_summary.append({
            'organism': organism,
            'protein_preserved': validation['proteins_match'],
            'nucleotide_identity': validation['nucleotide_identity'],
            'codons_changed': validation['codons_changed'],
            'output_file': output_path
        })
    
    # Print summary
    print("\n" + "="*80)
    print("OPTIMIZATION SUMMARY")
    print("="*80)
    print(f"\n{'Organism':20} | {'Protein OK':12} | {'Nuc. ID %':12} | {'Codons Changed':15}")
    print("-" * 70)
    
    for result in results_summary:
        status = "✅ YES" if result['protein_preserved'] else "❌ NO"
        print(f"{result['organism']:20} | {status:12} | {result['nucleotide_identity']:11.1f}% | {result['codons_changed']:15}")
    
    print("\n✅ All sequences optimized and saved to: data/modified_sequences/")
    print("\nNext steps:")
    print("1. Validate structures with ColabFold")
    print("2. Test against BLAST screening with: python code/test_screening.py")
    print("="*80 + "\n")