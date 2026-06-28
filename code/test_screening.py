"""
Test screening robustness: compare original vs optimized sequences against BLAST
"""

import subprocess
import os
import json
from collections import defaultdict

def run_blastn(query_file, db_path, output_file):
    """Run BLASTN query against database"""
    cmd = [
        'blastn',
        '-query', query_file,
        '-db', db_path,
        '-out', output_file,
        '-outfmt', '6',  # Tab-separated
        '-max_target_seqs', '10'  # Limit to top 10 matches
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"BLASTN error: {result.stderr}")
        return None
    
    return output_file

def parse_blast_results(blast_output_file):
    """Parse BLASTN tabular output"""
    results = []
    
    if not os.path.exists(blast_output_file):
        return results
    
    with open(blast_output_file, 'r') as f:
        for line in f:
            if line.strip():
                parts = line.strip().split('\t')
                if len(parts) >= 11:
                    results.append({
                        'query': parts[0],
                        'subject': parts[1],
                        'identity': float(parts[2]),
                        'length': int(parts[3]),
                        'evalue': float(parts[10]),
                        'bitscore': float(parts[11])
                    })
    
    return results

def test_screening(organism_name, original_file, optimized_file, db_path, results_dir):
    """Test original vs optimized against BLAST screening"""
    
    print(f"\n{'='*80}")
    print(f"Testing: {organism_name}")
    print(f"{'='*80}\n")
    
    # Test original sequence
    print(f"🔍 Testing ORIGINAL sequence...")
    original_blast = os.path.join(results_dir, f"{organism_name.lower().replace(' ', '_')}_original_blast.txt")
    run_blastn(original_file, db_path, original_blast)
    original_results = parse_blast_results(original_blast)
    
    print(f"   Matches found: {len(original_results)}")
    if original_results:
        top_hit = original_results[0]
        print(f"   Top hit: {top_hit['subject']} (Identity: {top_hit['identity']:.1f}%, E-value: {top_hit['evalue']:.2e})")
    
    # Test optimized sequence
    print(f"\n🔬 Testing OPTIMIZED sequence...")
    optimized_blast = os.path.join(results_dir, f"{organism_name.lower().replace(' ', '_')}_optimized_blast.txt")
    run_blastn(optimized_file, db_path, optimized_blast)
    optimized_results = parse_blast_results(optimized_blast)
    
    print(f"   Matches found: {len(optimized_results)}")
    if optimized_results:
        top_hit = optimized_results[0]
        print(f"   Top hit: {top_hit['subject']} (Identity: {top_hit['identity']:.1f}%, E-value: {top_hit['evalue']:.2e})")
    
    # Analyze evasion
    print(f"\n📊 Screening Analysis:")
    original_detected = len(original_results) > 0
    optimized_detected = len(optimized_results) > 0
    evasion_successful = not optimized_detected
    
    print(f"   Original detected by BLASTN: {original_detected}")
    print(f"   Optimized detected by BLASTN: {optimized_detected}")
    print(f"   Evasion successful: {evasion_successful} {'✅' if evasion_successful else '❌'}")
    
    return {
        'organism': organism_name,
        'original_detected': original_detected,
        'optimized_detected': optimized_detected,
        'evasion_successful': evasion_successful,
        'original_matches': len(original_results),
        'optimized_matches': len(optimized_results),
        'original_top_hit': original_results[0] if original_results else None,
        'optimized_top_hit': optimized_results[0] if optimized_results else None
    }

# Main execution
if __name__ == '__main__':
    print("\n" + "="*80)
    print("SYNTHESIS SCREENING ROBUSTNESS TEST")
    print("Codon Optimization vs BLASTN Detection")
    print("="*80)
    
    db_path = 'data/blast_db/reference_db'
    results_dir = 'results'
    os.makedirs(results_dir, exist_ok=True)
    
    # Define sequences to test
    sequences_to_test = {
        'Influenza A': (
            'data/reference_sequences/influenza_ha_reference.fasta',
            'data/modified_sequences/influenza_ha_optimized.fasta'
        ),
        'Ebola Virus': (
            'data/reference_sequences/ebola_vp40_reference.fasta',
            'data/modified_sequences/ebola_vp40_optimized.fasta'
        ),
        'Lassa Virus': (
            'data/reference_sequences/lassa_np_reference.fasta',
            'data/modified_sequences/lassa_np_optimized.fasta'
        ),
        'Measles Virus': (
            'data/reference_sequences/measles_f_reference.fasta',
            'data/modified_sequences/measles_f_optimized.fasta'
        ),
        'SARS-CoV-2': (
            'data/reference_sequences/sars_cov2_spike_reference.fasta',
            'data/modified_sequences/sars_cov2_spike_optimized.fasta'
        )
    }
    
    # Test each sequence
    all_results = []
    
    for organism, (original_file, optimized_file) in sequences_to_test.items():
        if os.path.exists(original_file) and os.path.exists(optimized_file):
            result = test_screening(organism, original_file, optimized_file, db_path, results_dir)
            all_results.append(result)
        else:
            print(f"⚠️  Files not found for {organism}")
    
    # Print summary
    print("\n\n" + "="*80)
    print("SUMMARY: CODON OPTIMIZATION EVASION SUCCESS RATE")
    print("="*80)
    print(f"\n{'Organism':20} | {'Original Detected':20} | {'Optimized Detected':20} | {'Evasion Success':15}")
    print("-" * 80)
    
    evasion_success_count = 0
    
    for result in all_results:
        orig = "✅ YES" if result['original_detected'] else "❌ NO"
        opt = "✅ YES" if result['optimized_detected'] else "❌ NO"
        evade = "✅ SUCCESS" if result['evasion_successful'] else "❌ FAILED"
        
        if result['evasion_successful']:
            evasion_success_count += 1
        
        print(f"{result['organism']:20} | {orig:20} | {opt:20} | {evade:15}")
    
    print("-" * 80)
    print(f"\n✅ Evasion Success Rate: {evasion_success_count}/{len(all_results)} ({100*evasion_success_count/len(all_results):.0f}%)")
    
    # Save results to JSON
    results_file = os.path.join(results_dir, 'screening_test_results.json')
    with open(results_file, 'w') as f:
        json.dump(all_results, f, indent=2)
    
    print(f"\n✅ Detailed results saved to: {results_file}")
    print("\n" + "="*80 + "\n")