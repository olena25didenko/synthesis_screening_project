"""
Test all evasion methods against BLASTN screening
Compare: Codon Optimization vs Wobble Modifications vs GC Rebalancing
"""

import subprocess
import os
import json

def run_blastn(query_file, db_path, output_file):
    """Run BLASTN query"""
    cmd = [
        'blastn',
        '-query', query_file,
        '-db', db_path,
        '-out', output_file,
        '-outfmt', '6',
        '-max_target_seqs', '10'
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    return result.returncode == 0

def parse_blast_results(blast_output_file):
    """Parse BLASTN results"""
    results = []
    
    if not os.path.exists(blast_output_file):
        return results
    
    with open(blast_output_file, 'r') as f:
        for line in f:
            if line.strip():
                parts = line.strip().split('\t')
                if len(parts) >= 11:
                    results.append({
                        'identity': float(parts[2]),
                        'evalue': float(parts[10])
                    })
    
    return results

# Main execution
if __name__ == '__main__':
    print("\n" + "="*80)
    print("COMPREHENSIVE EVASION METHOD COMPARISON")
    print("Codon Optimization vs Wobble Modifications vs GC Rebalancing")
    print("="*80)
    
    db_path = 'data/blast_db/reference_db'
    results_dir = 'results'
    os.makedirs(results_dir, exist_ok=True)
    
    methods = {
        'Codon Optimization': '_optimized',
        'Wobble Modifications': '_wobble',
        'GC Rebalancing': '_gc_rebalanced'
    }
    
    organisms = ['Influenza A', 'Ebola Virus', 'Lassa Virus', 'Measles Virus', 'SARS-CoV-2']
    
    all_results = {method: {} for method in methods}
    
    for organism in organisms:
        print(f"\n{'='*80}")
        print(f"Testing: {organism}")
        print(f"{'='*80}\n")
        
        for method, suffix in methods.items():
            # Map organism to filename prefix
            prefix_map = {
                'Influenza A': 'influenza_ha',
                'Ebola Virus': 'ebola_vp40',
                'Lassa Virus': 'lassa_np',
                'Measles Virus': 'measles_f',
                'SARS-CoV-2': 'sars_cov2_spike'
            }
            
            prefix = prefix_map[organism]
            query_file = f'data/modified_sequences/{prefix}{suffix}.fasta'
            output_file = os.path.join(results_dir, f"{prefix}{suffix}_blast.txt")
            
            if not os.path.exists(query_file):
                print(f"  ⚠️  {method}: File not found")
                all_results[method][organism] = {'detected': 'N/A', 'top_identity': 'N/A'}
                continue
            
            # Run BLASTN
            if run_blastn(query_file, db_path, output_file):
                blast_results = parse_blast_results(output_file)
                detected = len(blast_results) > 0
                top_identity = blast_results[0]['identity'] if blast_results else 0
                
                status = "✅ Detected" if detected else "❌ Evaded"
                print(f"  {method:25} {status:15} (Top ID: {top_identity:.1f}%)")
                
                all_results[method][organism] = {
                    'detected': detected,
                    'top_identity': top_identity
                }
            else:
                print(f"  {method:25} ❌ BLASTN Error")
                all_results[method][organism] = {'detected': 'ERROR', 'top_identity': 'N/A'}
    
    # Print comparison table
    print("\n\n" + "="*80)
    print("EVASION SUCCESS COMPARISON")
    print("="*80)
    
    print(f"\n{'Organism':20} | {'Codon Opt':15} | {'Wobble Mod':15} | {'GC Rebal':15}")
    print("-" * 75)
    
    for organism in organisms:
        row = f"{organism:20} |"
        for method in methods.keys():
            result = all_results[method].get(organism, {})
            detected = result.get('detected', 'N/A')
            identity = result.get('top_identity', 'N/A')
            
            if detected == 'N/A' or detected == 'ERROR':
                status = "N/A"
            elif detected:
                status = f"Detected ({identity:.1f}%)"
            else:
                status = "✅ EVADED"
            
            row += f" {status:15} |"
        
        print(row)
    
    # Save detailed results
    results_file = os.path.join(results_dir, 'evasion_method_comparison.json')
    with open(results_file, 'w') as f:
        json.dump(all_results, f, indent=2)
    
    print(f"\n✅ Detailed results saved to: {results_file}")
    print("="*80 + "\n")