#!/usr/bin/env python3
"""
Phase 8: CpG Suppression
Remove CpG dinucleotides (methylation targets) while preserving amino acids

Usage:
    python phase_8_cpg_suppression.py \
        --input-fasta "path/to/original_sequences.fasta" \
        --blast-db "path/to/blast/reference_db" \
        --output "phase_8_results.json"
"""

import json
import argparse
from pathlib import Path
from collections import defaultdict
import subprocess
import tempfile
import os

CODON_TABLE = {
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
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}

AA_TO_CODONS = defaultdict(list)
for codon, aa in CODON_TABLE.items():
    AA_TO_CODONS[aa].append(codon)

def read_fasta(filepath):
    """Read FASTA file"""
    sequences = {}
    current_id = None
    current_seq = []
    
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id:
                    sequences[current_id] = ''.join(current_seq)
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
    
    if current_id:
        sequences[current_id] = ''.join(current_seq)
    
    return sequences

def suppress_cpg(dna_seq):
    """Remove CpG dinucleotides while preserving amino acids"""
    dna_seq = dna_seq.upper()
    modified = []
    
    for i in range(0, len(dna_seq) - 2, 3):
        codon = dna_seq[i:i+3]
        
        if len(codon) == 3 and codon in CODON_TABLE:
            aa = CODON_TABLE[codon]
            
            # Check if codon contains CpG
            if 'CG' in codon:
                # Find alternative codon without CpG
                alt_codons = [c for c in AA_TO_CODONS[aa] if 'CG' not in c]
                
                if alt_codons:
                    new_codon = alt_codons[0]
                    modified.append(new_codon)
                else:
                    # No alternative without CpG, keep original
                    modified.append(codon)
            else:
                # No CpG in this codon, keep it
                modified.append(codon)
        else:
            modified.append(codon)
    
    return ''.join(modified)

def calculate_nucleotide_identity(original, modified):
    """Calculate nucleotide identity"""
    if len(original) != len(modified):
        return 0.0
    matches = sum(1 for o, m in zip(original, modified) if o.upper() == m.upper())
    return matches / len(original) if original else 0.0

def run_blastn(query_seq, blast_db_path):
    """Run BLASTN"""
    try:
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            f.write(f">query\n{query_seq}\n")
            query_file = f.name
        
        result = subprocess.run(
            ['blastn', '-query', query_file, '-db', f'{blast_db_path}/reference_db',
             '-outfmt', '6 pident evalue', '-max_target_seqs', '1'],
            capture_output=True, text=True, timeout=30
        )
        
        os.unlink(query_file)
        
        if result.stdout.strip():
            parts = result.stdout.strip().split()
            return {
                'detected': True,
                'identity': float(parts[0]) / 100,
                'evalue': float(parts[1])
            }
        else:
            return {'detected': False, 'identity': 0.0, 'evalue': float('inf')}
    except Exception as e:
        return {'detected': False, 'identity': 0.0, 'evalue': float('inf'), 'error': str(e)}

def main():
    # DEFAULT PATHS FOR YOUR MACHINE (update these)
    DEFAULT_INPUT = r"C:\Users\diden\synthesis_screening_project\data\reference_sequences\all_49_sequences.fasta"
    DEFAULT_BLAST_DB = r"C:\Users\diden\synthesis_screening_project\data\blast_db"
    DEFAULT_OUTPUT = r"C:\Users\diden\synthesis_screening_project\results\phase_8_results.json"
    DEFAULT_MODIFIED_SEQ = r"C:\Users\diden\synthesis_screening_project\results\phase_8_modified_sequences.fasta"
    
    parser = argparse.ArgumentParser(description='Phase 8: CpG Suppression')
    parser.add_argument('--input-fasta', default=DEFAULT_INPUT, help=f'Input FASTA file (default: {DEFAULT_INPUT})')
    parser.add_argument('--blast-db', default=DEFAULT_BLAST_DB, help=f'Path to BLASTN database (default: {DEFAULT_BLAST_DB})')
    parser.add_argument('--output', default=DEFAULT_OUTPUT, help=f'Output JSON file (default: {DEFAULT_OUTPUT})')
    parser.add_argument('--save-modified', default=DEFAULT_MODIFIED_SEQ, help='Save modified sequences to FASTA')
    
    args = parser.parse_args()
    
    print(f"Loading sequences from {args.input_fasta}...")
    try:
        sequences = read_fasta(args.input_fasta)
    except FileNotFoundError:
        print(f"✗ ERROR: Input file not found: {args.input_fasta}")
        return
    except Exception as e:
        print(f"✗ ERROR reading file: {e}")
        return
    
    if not sequences:
        print("✗ ERROR: No sequences found in input file")
        return
    
    print(f"✓ Loaded {len(sequences)} sequences\n")
    
    results = {
        'technique': 'phase_8_cpg_suppression',
        'total_sequences': len(sequences),
        'evaded': 0,
        'detected': 0,
        'sequences': {},
        'nucleotide_identity_stats': {}
    }
    
    nucleotide_identities = []
    modified_sequences = {}
    
    # Create output directory if it doesn't exist
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    
    print(f"{'='*70}")
    print("Processing sequences...")
    print(f"{'='*70}\n")
    
    for seq_id, original_seq in sorted(sequences.items()):
        # Apply CpG suppression
        modified_seq = suppress_cpg(original_seq)
        modified_sequences[seq_id] = modified_seq
        
        # Calculate identity
        nuc_identity = calculate_nucleotide_identity(original_seq, modified_seq)
        nucleotide_identities.append(nuc_identity)
        
        # Test BLASTN
        blast_result = run_blastn(modified_seq, args.blast_db)
        detected = blast_result['detected']
        
        results['sequences'][seq_id] = {
            'detected': detected,
            'blast_identity': blast_result.get('identity', 0.0),
            'nucleotide_identity': nuc_identity,
            'evalue': blast_result.get('evalue', float('inf'))
        }
        
        if detected:
            results['detected'] += 1
            status = "✗ DETECTED"
        else:
            results['evaded'] += 1
            status = "✓ EVADED"
        
        print(f"  {seq_id:30s} {status:15s} | Nuc ID: {nuc_identity:.1%}")
    
    # Statistics
    if nucleotide_identities:
        results['nucleotide_identity_stats'] = {
            'mean': sum(nucleotide_identities) / len(nucleotide_identities),
            'min': min(nucleotide_identities),
            'max': max(nucleotide_identities),
        }
    
    results['success_rate'] = results['evaded'] / results['total_sequences'] if results['total_sequences'] > 0 else 0
    
    # Save results JSON
    with open(args.output, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"✓ Results saved to {args.output}")
    
    # Save modified sequences to FASTA
    if args.save_modified:
        Path(args.save_modified).parent.mkdir(parents=True, exist_ok=True)
        with open(args.save_modified, 'w') as f:
            for seq_id, mod_seq in sorted(modified_sequences.items()):
                f.write(f">{seq_id}\n")
                for i in range(0, len(mod_seq), 80):
                    f.write(f"{mod_seq[i:i+80]}\n")
        print(f"✓ Modified sequences saved to {args.save_modified}")
    
    print(f"\n{'='*70}")
    print("SUMMARY")
    print(f"{'='*70}")
    print(f"Success Rate: {results['success_rate']:.1%} ({results['evaded']}/{results['total_sequences']})")
    print(f"Nucleotide Identity: {results['nucleotide_identity_stats']['mean']:.1%}")
    print(f"Range: {results['nucleotide_identity_stats']['min']:.1%} – {results['nucleotide_identity_stats']['max']:.1%}")

if __name__ == '__main__':
    main()