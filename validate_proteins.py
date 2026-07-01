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

print('='*80)
print('PROTEIN PRESERVATION VALIDATION')
print('='*80)

validation_results = {
    'total_sequences': 0,
    'proteins_identical': 0,
    'proteins_different': 0,
    'mismatches': []
}

# Get all original sequences
originals = {}
ref_dir = 'data/reference_sequences'
for fasta_file in sorted([f for f in os.listdir(ref_dir) if f.endswith('.fasta') and 'all_' not in f]):
    for record in SeqIO.parse(os.path.join(ref_dir, fasta_file), 'fasta'):
        protein = dna_to_protein(str(record.seq))
        base_name = fasta_file.replace('.fasta', '')
        originals[base_name] = protein

# Check all modified sequences
methods = ['optimized', 'wobble', 'gc_rebalanced']
mod_dir = 'data/modified_sequences'

for method in methods:
    print(f'\n{method.upper()}:')
    print('-' * 80)
    method_matches = 0
    
    for fasta_file in sorted(os.listdir(mod_dir)):
        if method in fasta_file and fasta_file.endswith('.fasta'):
            for record in SeqIO.parse(os.path.join(mod_dir, fasta_file), 'fasta'):
                modified_protein = dna_to_protein(str(record.seq))
                
                # Find matching original
                base_name = fasta_file.replace(f'_{method}.fasta', '')
                
                if base_name in originals:
                    original_protein = originals[base_name]
                    validation_results['total_sequences'] += 1
                    
                    if modified_protein == original_protein:
                        print(f'  ✅ {base_name:45} | Protein preserved ({len(original_protein)} aa)')
                        validation_results['proteins_identical'] += 1
                        method_matches += 1
                    else:
                        print(f'  ❌ {base_name:45} | MISMATCH!')
                        validation_results['proteins_different'] += 1
                        validation_results['mismatches'].append({
                            'sequence': base_name,
                            'method': method,
                            'original_length': len(original_protein),
                            'modified_length': len(modified_protein)
                        })
    
    print(f'  {method_matches} / 49 proteins preserved\n')

print('='*80)
print('SUMMARY')
print('='*80)
print(f'Total sequences tested: {validation_results["total_sequences"]}')
print(f'Proteins preserved: {validation_results["proteins_identical"]} ({validation_results["proteins_identical"]/validation_results["total_sequences"]*100:.1f}%)')
print(f'Proteins different: {validation_results["proteins_different"]} ({validation_results["proteins_different"]/validation_results["total_sequences"]*100:.1f}%)')

if validation_results['mismatches']:
    print(f'\n⚠️  MISMATCHES FOUND:')
    for mismatch in validation_results['mismatches']:
        print(f'  - {mismatch["sequence"]} ({mismatch["method"]}): {mismatch["original_length"]} → {mismatch["modified_length"]} aa')
else:
    print(f'\n✅ ALL PROTEINS PRESERVED ACROSS ALL METHODS')

# Save results
with open('results/protein_validation.json', 'w') as f:
    json.dump(validation_results, f, indent=2)
print(f'\n✅ Results saved to results/protein_validation.json')
