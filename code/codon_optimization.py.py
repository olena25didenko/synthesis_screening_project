"""
Codon optimization - AUTO-DISCOVERS all FASTA files in reference_sequences/
"""

from Bio import SeqIO
import json
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

CODONS_FOR_AA = {}
for codon, aa in GENETIC_CODE.items():
    if aa not in CODONS_FOR_AA:
        CODONS_FOR_AA[aa] = []
    CODONS_FOR_AA[aa].append(codon)

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

def optimize_sequence(dna_seq, optimization_map):
    optimized = []
    for i in range(0, len(dna_seq) - 2, 3):
        codon = dna_seq[i:i+3].upper()
        if codon in optimization_map:
            optimized.append(optimization_map[codon])
        else:
            optimized.append(codon)
    return ''.join(optimized)

def validate_optimization(original_seq, optimized_seq):
    original_protein = dna_to_protein(original_seq)
    optimized_protein = dna_to_protein(optimized_seq)
    
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

if __name__ == '__main__':
    print("\n" + "="*80)
    print("CODON OPTIMIZATION FOR SYNTHESIS SCREENING PROJECT")
    print("AUTO-DISCOVERING ALL SEQUENCES")
    print("="*80 + "\n")
    
    try:
        with open('data/codon_usage/optimization_maps.json', 'r') as f:
            all_maps = json.load(f)
    except FileNotFoundError:
        print("❌ optimization_maps.json not found. Run extract_codon_usage.py first!")
        exit(1)
    
    os.makedirs('data/modified_sequences', exist_ok=True)
    
    # AUTO-DISCOVER all FASTA files
    ref_dir = 'data/reference_sequences'
    # NEW (skips combined file, only processes individual references):
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
        
        # Find best matching organism in optimization maps
        best_organism = None
        for organism in all_maps.keys():
            if organism.lower().replace(' ', '_') in fasta_file.lower():
                best_organism = organism
                break
        
        if not best_organism:
            # Fallback: use first organism if can't match
            best_organism = list(all_maps.keys())[0]
        
        optimization_map = all_maps[best_organism]
        optimized_seq = optimize_sequence(original_seq, optimization_map)
        validation = validate_optimization(original_seq, optimized_seq)
        
        # Save
        output_filename = fasta_file.replace('.fasta', '_optimized.fasta')
        output_path = os.path.join('data/modified_sequences', output_filename)
        
        with open(output_path, 'w') as f:
            f.write(f">{record.id}_codon_optimized\n")
            f.write(optimized_seq + "\n")
        
        print(f"✅ ({validation['nucleotide_identity']:.1f}% identity)")
        
        results_summary.append({
            'file': fasta_file,
            'protein_preserved': validation['proteins_match'],
            'nucleotide_identity': validation['nucleotide_identity'],
            'codons_changed': validation['codons_changed'],
            'output_file': output_path
        })
    
    # Print summary
    print("\n" + "="*80)
    print(f"CODON OPTIMIZATION SUMMARY ({len(results_summary)} sequences)")
    print("="*80)
    print(f"\n{'File':40} | {'Protein OK':12} | {'Nuc. ID %':12}")
    print("-" * 70)
    
    for result in results_summary:
        status = "✅ YES" if result['protein_preserved'] else "❌ NO"
        print(f"{result['file']:40} | {status:12} | {result['nucleotide_identity']:11.1f}%")
    
    print(f"\n✅ Processed {len(results_summary)} sequences")
    print("="*80 + "\n")