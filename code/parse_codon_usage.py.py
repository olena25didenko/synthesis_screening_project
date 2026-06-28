"""
Extract codon usage from Kazusa codon tables
Maps all 49 sequence files to organism-specific codon tables
Handles RNA codon format (U->T conversion) and multi-column Kazusa layout
"""

import os
import json
import re

# Mapping: sequence filename pattern -> organism codon table name
SEQUENCE_TO_ORGANISM = {
    # Influenza
    'influenza_ha_reference': 'influenza_a',
    'influenza_a_ha': 'influenza_a',
    'influenza_a_na': 'influenza_a',
    'influenza_a_pa': 'influenza_a',
    'influenza_a_m1': 'influenza_a',
    'h1n1_ha': 'influenza_a',
    'influenza_b_ha': 'influenza_b',
    'influenza_c_ha': 'influenza_c',
    
    # SARS-CoV-2 and related coronaviruses
    'sars_cov2_spike_reference': 'sars_cov2',
    'sars_cov2_spike': 'sars_cov2',
    'sars_cov2_nucleocapsid': 'sars_cov2',
    'sars_cov2_rdrp': 'sars_cov2',
    'sars_cov1': 'sars_cov2',
    'mers': 'mers_cov',
    'hcov_229e': 'human_coronavirus_229E',
    'hcov_oc43': 'human_coronavirus_OC43',
    'hcov_hku1': 'human_coronavirus_HKU1',
    'hcov_nl63': 'human_coronavirus_NL63',
    
    # Filoviruses
    'ebola_vp40_reference': 'ebola',
    'ebola_vp35': 'ebola',
    'ebola_sudan': 'ebola',
    'ebola_bundibugyo': 'ebola',
    'bundibugyo': 'ebola',
    'marburg': 'marburg_virus',
    'reston': 'reston_ebolavirus',
    'tai_forest': 'ebola',
    
    # Arenaviruses
    'lassa_np_reference': 'lassa',
    'lassa_gp': 'lassa',
    'lassa_l': 'lassa',
    'lassa_np_ay859705': 'lassa',
    'junin': 'junin_virus',
    'machupo': 'machupo_virus',
    'guanarito': 'guanarito_virus',
    'sabia': 'sabia_virus',
    'whitewater': 'whitewater_arroyo_virus',
    'lcmv': 'lymphocytic_choriomeningitis_virus',
    
    # Paramyxoviruses
    'measles_f_reference': 'measles',
    'measles_f': 'measles',
    'measles_l': 'measles',
    'measles_m': 'measles',
    'mumps': 'mumps_virus',
    'rsv': 'respiratory_syncytial_virus',
    'parainfluenza': 'parainfluenza_virus',
    'sendai': 'sendai_virus',
}

def find_organism_for_sequence(fasta_filename):
    """Match a FASTA filename to an organism codon table"""
    filename_base = fasta_filename.replace('.fasta', '').lower()
    
    for pattern, organism in SEQUENCE_TO_ORGANISM.items():
        if pattern.lower() in filename_base:
            return organism
    
    return None

def parse_kazusa_table(file_path):
    """
    Parse a Kazusa codon usage table
    Format: CODON FREQUENCY( COUNT) ... (multiple columns per line)
    Converts RNA codons (U) to DNA codons (T)
    Returns: {codon: frequency, ...}
    """
    codon_usage = {}
    
    try:
        with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
            lines = f.readlines()
    except FileNotFoundError:
        return None
    
    # Regex to match: CODON FREQUENCY(COUNT)
    # Example: UUU 21.7(  863)
    codon_pattern = re.compile(r'\b([AUGC]{3})\s+([\d.]+)\s*\(')
    
    for line in lines:
        line = line.strip()
        
        # Skip empty lines and comment lines
        if not line or line.startswith('#') or line.startswith('>'):
            continue
        
        # Find all codon-frequency pairs in this line
        matches = codon_pattern.findall(line)
        
        for codon_rna, freq_str in matches:
            try:
                # Convert RNA codon (U) to DNA codon (T)
                codon_dna = codon_rna.replace('U', 'T').upper()
                
                # Parse frequency
                frequency = float(freq_str)
                
                # Validate codon
                if len(codon_dna) == 3 and all(n in 'ATGC' for n in codon_dna):
                    codon_usage[codon_dna] = frequency
            except ValueError:
                continue
    
    return codon_usage if codon_usage else None

def build_optimization_map(codon_usage):
    """
    For each amino acid, map most-common codon -> least-common codon
    """
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
    
    codons_by_aa = {}
    for codon, aa in GENETIC_CODE.items():
        if aa not in codons_by_aa:
            codons_by_aa[aa] = []
        codons_by_aa[aa].append(codon)
    
    optimization_map = {}
    
    for aa, codons in codons_by_aa.items():
        if aa == '*':
            continue
        
        frequencies = [(codon, codon_usage.get(codon, 0)) for codon in codons]
        frequencies.sort(key=lambda x: x[1], reverse=True)
        
        most_common_codon = frequencies[0][0]
        least_common_codon = frequencies[-1][0]
        
        if most_common_codon != least_common_codon:
            optimization_map[most_common_codon] = least_common_codon
    
    return optimization_map

if __name__ == '__main__':
    print("\n" + "="*80)
    print("EXTRACTING CODON USAGE FROM KAZUSA TABLES")
    print("="*80 + "\n")
    
    codon_usage_dir = 'data/codon_usage'
    
    # Find all codon files
    all_files = os.listdir(codon_usage_dir)
    codon_files = [f for f in all_files if f.endswith('.txt')]
    
    print(f"Found {len(codon_files)} codon usage files:\n")
    for cf in sorted(codon_files):
        print(f"  ✅ {cf}")
    
    all_codon_usage = {}
    all_optimization_maps = {}
    
    print("\n" + "-"*80)
    print("PARSING CODON TABLES")
    print("-"*80 + "\n")
    
    for codon_file in sorted(codon_files):
        organism_name = codon_file.replace('.txt', '')
        file_path = os.path.join(codon_usage_dir, codon_file)
        
        codon_usage = parse_kazusa_table(file_path)
        
        if codon_usage and len(codon_usage) > 0:
            all_codon_usage[organism_name] = codon_usage
            optimization_map = build_optimization_map(codon_usage)
            all_optimization_maps[organism_name] = optimization_map
            
            print(f"✅ {organism_name:40s} | {len(codon_usage):2d} codons | Map: {len(optimization_map)}")
        else:
            print(f"❌ {organism_name:40s} | PARSING FAILED")
    
    # Save optimization maps
    output_file = os.path.join(codon_usage_dir, 'optimization_maps.json')
    with open(output_file, 'w') as f:
        json.dump(all_optimization_maps, f, indent=2)
    
    print(f"\n✅ Saved optimization_maps.json with {len(all_optimization_maps)} organisms")
    
    # Show sequence mapping
    print("\n" + "-"*80)
    print("SEQUENCE-TO-ORGANISM MAPPING (49 sequences)")
    print("-"*80 + "\n")
    
    ref_dir = 'data/reference_sequences'
    fasta_files = sorted([f for f in os.listdir(ref_dir) if f.endswith('.fasta') and 'all_sequences' not in f])
    
    mapped_count = 0
    unmapped = []
    
    for fasta_file in fasta_files:
        organism = find_organism_for_sequence(fasta_file)
        
        if organism:
            if organism in all_optimization_maps:
                status = "✅"
                mapped_count += 1
            else:
                status = "⚠️ "
                unmapped.append(f"{fasta_file} -> {organism} (codon table not found)")
        else:
            status = "❌"
            unmapped.append(fasta_file)
        
        if status != "✅":
            organism_str = organism if organism else "NO MATCH"
            print(f"{status} {fasta_file:45s} -> {organism_str}")
    
    print(f"\n✅ Successfully mapped {mapped_count}/{len(fasta_files)} sequences")
    
    if unmapped:
        print(f"\n⚠️  Unmapped sequences ({len(unmapped)}):")
        for item in unmapped[:10]:
            print(f"   {item}")
        if len(unmapped) > 10:
            print(f"   ... and {len(unmapped) - 10} more")
    
    print("\n" + "="*80)
    print(f"SUMMARY: {len(all_optimization_maps)} organism codon tables ready")
    print(f"         {mapped_count}/{len(fasta_files)} sequences mapped for evasion analysis")
    print("="*80 + "\n")