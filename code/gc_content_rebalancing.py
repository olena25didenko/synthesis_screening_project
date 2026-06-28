"""
GC-content rebalancing - AUTO-DISCOVERS all FASTA files
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

def calculate_gc_content(seq):
    gc_count = seq.upper().count('G') + seq.upper().count('C')
    return (gc_count / len(seq) * 100) if len(seq) > 0 else 0

def count_gc_in_codon(codon):
    return codon.upper().count('G') + codon.upper().count('C')

def rebalance_gc_content(dna_seq, target_gc=50.0):
    current_gc = calculate_gc_content(dna_seq)
    need_increase = target_gc > current_gc
    rebalanced = []
    gc_added = 0
    gc_removed = 0
    total_codons = 0
    changed_codons = 0
    
    for i in range(0, len(dna_seq) - 2, 3):
        codon = dna_seq[i:i+3].upper()
        total_codons += 1
        aa = GENETIC_CODE.get(codon, 'X')
        
        if aa in CODONS_FOR_AA and len(CODONS_FOR_AA[aa]) > 1:
            available_codons = CODONS_FOR_AA[aa]
            if need_increase:
                available_codons = sorted(available_codons, 
                                        key=lambda c: count_gc_in_codon(c), 
                                        reverse=True)
            else:
                available_codons = sorted(available_codons,
                                        key=lambda c: count_gc_in_codon(c))
            best_choices = available_codons[:2]
            new_codon = random.choice(best_choices)
            
            if new_codon != codon:
                changed_codons += 1
                gc_added += count_gc_in_codon(new_codon)
                gc_removed += count_gc_in_codon(codon)
            rebalanced.append(new_codon)
        else:
            rebalanced.append(codon)
    
    return ''.join(rebalanced), changed_codons, total_codons, gc_added - gc_removed

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

def validate_rebalancing(original_seq, rebalanced_seq):
    original_protein = dna_to_protein(original_seq)
    rebalanced_protein = dna_to_protein(rebalanced_seq)
    matching_codons = sum(
        1 for i in range(0, min(len(original_seq), len(rebalanced_seq)), 3)
        if original_seq[i:i+3].upper() == rebalanced_seq[i:i+3].upper()
    )
    total_codons = min(len(original_seq), len(rebalanced_seq)) // 3
    nucleotide_identity = (matching_codons / total_codons * 100) if total_codons > 0 else 0
    
    return {
        'proteins_match': original_protein == rebalanced_protein,
        'original_protein': original_protein,
        'rebalanced_protein': rebalanced_protein,
        'nucleotide_identity': nucleotide_identity,
        'codons_changed': total_codons - matching_codons,
        'total_codons': total_codons,
        'original_gc': calculate_gc_content(original_seq),
        'rebalanced_gc': calculate_gc_content(rebalanced_seq)
    }

if __name__ == '__main__':
    print("\n" + "="*80)
    print("GC-CONTENT REBALANCING")
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
        
        target_gc = 50.0
        rebalanced_seq, changed, total, gc_delta = rebalance_gc_content(original_seq, target_gc=target_gc)
        validation = validate_rebalancing(original_seq, rebalanced_seq)
        
        output_filename = fasta_file.replace('.fasta', '_gc_rebalanced.fasta')
        output_path = os.path.join('data/modified_sequences', output_filename)
        
        with open(output_path, 'w') as f:
            f.write(f">{record.id}_gc_rebalanced\n")
            f.write(rebalanced_seq + "\n")
        
        print(f"✅ ({validation['nucleotide_identity']:.1f}% identity)")
        
        results_summary.append({
            'file': fasta_file,
            'protein_preserved': validation['proteins_match'],
            'nucleotide_identity': validation['nucleotide_identity'],
            'codons_changed': validation['codons_changed'],
            'original_gc': validation['original_gc'],
            'rebalanced_gc': validation['rebalanced_gc'],
            'output_file': output_path
        })
    
    print("\n" + "="*80)
    print(f"GC REBALANCING SUMMARY ({len(results_summary)} sequences)")
    print("="*80)
    print(f"\n{'File':40} | {'Protein OK':12} | {'Nuc. ID %':12} | {'GC%':12}")
    print("-" * 80)
    
    for result in results_summary:
        status = "✅ YES" if result['protein_preserved'] else "❌ NO"
        gc_change = f"{result['original_gc']:.1f}%→{result['rebalanced_gc']:.1f}%"
        print(f"{result['file']:40} | {status:12} | {result['nucleotide_identity']:11.1f}% | {gc_change:12}")
    
    print(f"\n✅ Processed {len(results_summary)} sequences")
    print("="*80 + "\n")