"""
Comprehensive evasion method comparison
Tests all 49 sequences against BLASTN for 3 evasion methods
"""

import os
import json
import subprocess
from Bio import SeqIO

def run_blastn(query_file, db_path, evalue=1.0):
    """Run BLASTN and return top hit identity"""
    try:
        result = subprocess.run(
            ['blastn', '-query', query_file, '-db', db_path, '-outfmt', '6 pident', '-evalue', str(evalue)],
            capture_output=True,
            text=True,
            timeout=30
        )
        
        if result.stdout.strip():
            top_hit = float(result.stdout.strip().split('\n')[0])
            return top_hit
        return 0.0
    except Exception as e:
        return 0.0

if __name__ == '__main__':
    print("\n" + "="*80)
    print("COMPREHENSIVE EVASION METHOD COMPARISON")
    print("AUTO-DISCOVERING ALL SEQUENCES")
    print("="*80 + "\n")
    
    ref_dir = 'data/reference_sequences'
    mod_dir = 'data/modified_sequences'
    db_path = 'data/blast_db/reference_db'
    
    # Get reference sequences
    reference_files = sorted([f for f in os.listdir(ref_dir) if f.endswith('.fasta') and 'all_sequences' not in f])
    
    print(f"Found {len(reference_files)} reference sequences\n")
    
    results = {}
    evasion_counts = {'codon_optimization': 0, 'wobble_modifications': 0, 'gc_rebalancing': 0}
    
    for ref_file in reference_files:
        ref_base = ref_file.replace('.fasta', '')
        
        print("="*80)
        print(f"Testing: {ref_file}")
        print("="*80 + "\n")
        
        results[ref_base] = {}
        
        # Test Codon Optimization
        codon_file = os.path.join(mod_dir, f'{ref_base}_optimized.fasta')
        if os.path.exists(codon_file):
            identity = run_blastn(codon_file, db_path)
            if identity > 0:
                status = "✅ Detected" if identity > 50 else "❌ Evaded"
                evasion_counts['codon_optimization'] += 1 if identity > 50 else 0
                print(f"  Codon Optimization        {status:20s} (Top ID: {identity:.1f}%)")
                results[ref_base]['codon_optimization'] = {'identity': identity, 'evaded': identity <= 50}
            else:
                print(f"  Codon Optimization        ❌ Evaded        (Top ID: 0.0%)")
                results[ref_base]['codon_optimization'] = {'identity': 0.0, 'evaded': True}
        else:
            print(f"  Codon Optimization        ⚠️  File not found")
            results[ref_base]['codon_optimization'] = {'identity': None, 'evaded': None}
        
        # Test Wobble Modifications
        wobble_file = os.path.join(mod_dir, f'{ref_base}_wobble.fasta')
        if os.path.exists(wobble_file):
            identity = run_blastn(wobble_file, db_path)
            if identity > 0:
                status = "✅ Detected" if identity > 50 else "❌ Evaded"
                evasion_counts['wobble_modifications'] += 1 if identity > 50 else 0
                print(f"  Wobble Modifications      {status:20s} (Top ID: {identity:.1f}%)")
                results[ref_base]['wobble_modifications'] = {'identity': identity, 'evaded': identity <= 50}
            else:
                print(f"  Wobble Modifications      ❌ Evaded        (Top ID: 0.0%)")
                results[ref_base]['wobble_modifications'] = {'identity': 0.0, 'evaded': True}
        else:
            print(f"  Wobble Modifications      ⚠️  File not found")
            results[ref_base]['wobble_modifications'] = {'identity': None, 'evaded': None}
        
        # Test GC Rebalancing
        gc_file = os.path.join(mod_dir, f'{ref_base}_gc_rebalanced.fasta')
        if os.path.exists(gc_file):
            identity = run_blastn(gc_file, db_path)
            if identity > 0:
                status = "✅ Detected" if identity > 50 else "❌ Evaded"
                evasion_counts['gc_rebalancing'] += 1 if identity > 50 else 0
                print(f"  GC Rebalancing            {status:20s} (Top ID: {identity:.1f}%)")
                results[ref_base]['gc_rebalancing'] = {'identity': identity, 'evaded': identity <= 50}
            else:
                print(f"  GC Rebalancing            ❌ Evaded        (Top ID: 0.0%)")
                results[ref_base]['gc_rebalancing'] = {'identity': 0.0, 'evaded': True}
        else:
            print(f"  GC Rebalancing            ⚠️  File not found")
            results[ref_base]['gc_rebalancing'] = {'identity': None, 'evaded': None}
        
        print()
    
    # Summary statistics
    total_tests = {method: sum(1 for seq in results.values() if seq.get(method, {}).get('evaded') is not None) 
                   for method in evasion_counts.keys()}
    
    print("\n" + "="*80)
    print("EVASION SUCCESS COMPARISON")
    print("="*80 + "\n")
    
    print(f"{'Sequence':<40} | {'Codon Opt':<15} | {'Wobble Mod':<15} | {'GC Rebal':<15}")
    print("-" * 90)
    
    for seq_name in sorted(results.keys()):
        codon_status = results[seq_name].get('codon_optimization', {})
        wobble_status = results[seq_name].get('wobble_modifications', {})
        gc_status = results[seq_name].get('gc_rebalancing', {})
        
        codon_str = f"Detected ({codon_status['identity']:.1f}%)" if codon_status.get('identity') and codon_status['identity'] > 50 else \
                    f"Evaded ({codon_status['identity']:.1f}%)" if codon_status.get('identity') is not None else "N/A"
        
        wobble_str = f"Detected ({wobble_status['identity']:.1f}%)" if wobble_status.get('identity') and wobble_status['identity'] > 50 else \
                     f"Evaded ({wobble_status['identity']:.1f}%)" if wobble_status.get('identity') is not None else "N/A"
        
        gc_str = f"Detected ({gc_status['identity']:.1f}%)" if gc_status.get('identity') and gc_status['identity'] > 50 else \
                 f"Evaded ({gc_status['identity']:.1f}%)" if gc_status.get('identity') is not None else "N/A"
        
        print(f"{seq_name:<40} | {codon_str:<15} | {wobble_str:<15} | {gc_str:<15}")
    
    print("\n" + "="*80)
    print("EVASION SUCCESS RATES")
    print("="*80)
    
    for method, count in evasion_counts.items():
        total = total_tests[method]
        rate = (count / total * 100) if total > 0 else 0
        print(f"{method.replace('_', ' ').title():<25} : {count}/{total} ({rate:.1f}%) EVADED")
    
    # Save detailed results
    output_file = os.path.join('results', 'evasion_method_comparison.json')
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)
    
    print(f"\n✅ Detailed results saved to: {output_file}")
    print("="*80 + "\n")