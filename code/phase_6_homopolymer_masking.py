#!/usr/bin/env python3
"""
Phase 6: Repeat & Homopolymer Masking
Break homopolymer runs (AAAA → AACA) to prevent BLASTN ~11bp exact-match anchoring

Usage:
    python phase_6_homopolymer_masking.py \
        --input-fasta "path/to/original_sequences.fasta" \
        --blast-db "path/to/blast/reference_db" \
        --output "phase_6_results.json"
"""

import json
import argparse
from pathlib import Path
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

AA_TO_CODONS = {
    'F': ['TTT', 'TTC'], 'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
    'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'], 'Y': ['TAT', 'TAC'],
    'C': ['TGT', 'TGC'], 'W': ['TGG'], 'P': ['CCT', 'CCC', 'CCA', 'CCG'],
    'H': ['CAT', 'CAC'], 'Q': ['CAA', 'CAG'], 'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
    'I': ['ATT', 'ATC', 'ATA'], 'M': ['ATG'], 'T': ['ACT', 'ACC', 'ACA', 'ACG'],
    'N': ['AAT', 'AAC'], 'K': ['AAA', 'AAG'], 'V': ['GTT', 'GTC', 'GTA', 'GTG'],
    'A': ['GCT', 'GCC', 'GCA', 'GCG'], 'D': ['GAT', 'GAC'], 'E': ['GAA', 'GAG'],
    'G': ['GGT', 'GGC', 'GGA', 'GGG'], '*': ['TAA', 'TAG', 'TGA']
}

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

def mask_homopolymers(dna_seq):
    """Break homopolymer runs to prevent BLASTN anchor-based detection"""
    modified = list(dna_seq.upper())
    i = 0
    
    while i < len(modified) - 10:
        # Check for homopolymer run (same base ≥10 times)
        if len(set(modified[i:i+10])) == 1:
            # Found homopolymer, insert alternative base at position i+9
            base = modified[i]
            alt_bases = [b for b in 'ATGC' if b != base]
            
            # Get codon position
            codon_pos = (i // 3) * 3
            in_codon_offset = (i + 9) % 3
            
            # Insert alternative base
            if alt_bases:
                modified.insert(i + 9, alt_bases[0])
            i += 10
        else:
            i += 1
    
    return ''.join(modified)

def calculate_nucleotide_identity(original, modified):
    """Calculate nucleotide identity"""
    if not original or not modified:
        return 0.0
    
    # Align sequences and calculate identity
    matches = 0
    min_len = min(len(original), len(modified))
    
    for i in range(min_len):
        if original[i].upper() == modified[i].upper():
            matches += 1
    
    return matches / max(len(original), len(modified))

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
    DEFAULT_OUTPUT = r"C:\Users\diden\synthesis_screening_project\results\phase_6_results.json"
    DEFAULT_MODIFIED_SEQ = r"C:\Users\diden\synthesis_screening_project\results\phase_6_modified_sequences.fasta"
    
    parser = argparse.ArgumentParser(description='Phase 6: Homopolymer Masking')
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
        'technique': 'phase_6_homopolymer_masking',
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
        # Apply masking
        modified_seq = mask_homopolymers(original_seq)
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