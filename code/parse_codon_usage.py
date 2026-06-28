"""
Extract codon usage statistics from Kazusa codon usage tables
Organize by amino acid and identify rare vs common codons for optimization
"""

from collections import defaultdict
import os
import json

# Standard genetic code (codon -> amino acid)
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

def parse_kazusa_line(line):
    """
    Parse a Kazusa codon usage line
    Format: "TTT  17.8(293396)  TTC  20.2(332322)  TTA  8.1(133967)  TTG  14.6(240030)"
    Returns: list of (codon, frequency) tuples
    """
    codons_data = []
    parts = line.split()
    
    i = 0
    while i < len(parts):
        if len(parts[i]) == 3 and parts[i].isupper():  # It's a codon
            codon = parts[i]
            if i + 1 < len(parts):
                # Extract frequency from "17.8(293396)" format
                freq_str = parts[i + 1]
                try:
                    # Get frequency per thousand (before parenthesis)
                    freq = float(freq_str.split('(')[0])
                    codons_data.append((codon, freq))
                except (ValueError, IndexError):
                    pass
            i += 2
        else:
            i += 1
    
    return codons_data

def extract_from_kazusa_file(filename):
    """
    Read Kazusa codon usage file and organize by amino acid
    """
    codon_freqs = defaultdict(list)
    
    try:
        with open(filename, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                
                # Skip headers and empty lines
                if not line or 'fields:' in line or 'Coding GC' in line or 'Format:' in line:
                    continue
                
                # Parse codon data
                codons_data = parse_kazusa_line(line)
                
                # Organize by amino acid
                for codon, freq in codons_data:
                    aa = GENETIC_CODE.get(codon, 'X')
                    codon_freqs[aa].append((codon, freq))
        
        # Sort each amino acid's codons by frequency
        for aa in codon_freqs:
            codon_freqs[aa].sort(key=lambda x: x[1], reverse=True)
        
        return codon_freqs
    
    except Exception as e:
        print(f"Error reading {filename}: {e}")
        return None

def print_codon_summary(organism_name, codon_freqs):
    """Print organized codon usage"""
    print(f"\n{'='*80}")
    print(f"Codon Usage: {organism_name.upper()}")
    print(f"{'='*80}\n")
    print(f"{'AA':3} | {'Most Common (Rare->Optimize)':30} | All Codons (freq/1000)")
    print("-" * 100)
    
    for aa in sorted(codon_freqs.keys()):
        if aa == '*':
            continue  # Skip stop codons
        
        codons_list = codon_freqs[aa]
        most_common = f"{codons_list[0][0]}({codons_list[0][1]:.1f})" if codons_list else "N/A"
        least_common = f"{codons_list[-1][0]}({codons_list[-1][1]:.1f})" if len(codons_list) > 1 else most_common
        
        all_codons = ', '.join([f"{c[0]}({c[1]:.1f})" for c in codons_list])
        
        target = f"{most_common} -> {least_common}"
        print(f"{aa:3} | {target:30} | {all_codons}")
    
    print("\n")

def get_optimization_map(codon_freqs):
    """
    Create mapping: common codon -> rare codon (for optimization)
    This helps with "simple" optimization strategy
    """
    optimization_map = {}
    
    for aa in codon_freqs:
        if aa == '*':
            continue
        
        codons_list = codon_freqs[aa]
        if len(codons_list) > 1:
            # Map most common -> least common (maximize change)
            most_common = codons_list[0][0]
            least_common = codons_list[-1][0]
            optimization_map[most_common] = least_common
    
    return optimization_map

def save_optimization_maps(all_results):
    """Save optimization maps to JSON for easy import"""
    output = {}
    
    for organism, codon_freqs in all_results.items():
        opt_map = get_optimization_map(codon_freqs)
        output[organism] = opt_map
    
    # Save as JSON
    output_file = 'data/codon_usage/optimization_maps.json'
    with open(output_file, 'w') as f:
        json.dump(output, f, indent=2)
    
    print(f"\n✅ Optimization maps saved to: {output_file}")
    return output_file

# Main execution
if __name__ == '__main__':
    codon_dir = 'data/codon_usage'
    
    print("\n" + "="*80)
    print("EXTRACTING CODON USAGE FROM KAZUSA FILES")
    print("="*80)
    
    # Find all codon usage files
    files_to_process = {
        'Influenza A': 'influenza_a.txt',
        'Ebola Virus': 'ebola.txt',
        'Lassa Virus': 'lassa.txt',
        'Measles Virus': 'measles.txt',
        'SARS-CoV-2': 'sars_cov2.txt'
    }
    
    all_results = {}
    
    for organism_name, filename in files_to_process.items():
        filepath = os.path.join(codon_dir, filename)
        
        if os.path.exists(filepath):
            print(f"\n📖 Processing: {filename}...")
            codon_freqs = extract_from_kazusa_file(filepath)
            
            if codon_freqs:
                all_results[organism_name] = codon_freqs
                print_codon_summary(organism_name, codon_freqs)
            else:
                print(f"❌ Failed to parse {filename}")
        else:
            print(f"⚠️  File not found: {filepath}")
    
    # Save optimization maps
    if all_results:
        save_optimization_maps(all_results)
        
        print("\n" + "="*80)
        print("SUMMARY")
        print("="*80)
        print(f"✅ Successfully processed {len(all_results)} organisms")
        print(f"✅ Optimization maps ready for codon_optimization.py")
        print("\nNext steps:")
        print("1. Run: python code/codon_optimization.py")
        print("2. Test optimized sequences against BLAST screening")
        print("="*80 + "\n")