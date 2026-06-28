"""
Download codon usage tables from Kazusa for all 49 organisms
"""

import os
from Bio import Entrez

# Set your email (NCBI requires this)
Entrez.email = "didenkoolena2504@gmail.com"

# Mapping: filename -> Kazusa search term
organisms_to_download = {
    # Influenza
    'influenza_a_ha_cy163680': 'Influenza A virus',
    'influenza_a_ha_jq708657': 'Influenza A virus',
    'influenza_a_ha_kf821935': 'Influenza A virus',
    'influenza_a_na_cy163681': 'Influenza A virus',
    'influenza_a_pa_cy163684': 'Influenza A virus',
    'influenza_a_m1_cy163686': 'Influenza A virus',
    'influenza_b_ha_lc107905': 'Influenza B virus',
    'influenza_c_ha_ab211639': 'Influenza C virus',
    'h1n1_ha_gq152857': 'Influenza A virus',
    
    # Coronaviruses
    'sars_cov2_spike_ov109917': 'SARS-CoV-2',
    'sars_cov2_nucleocapsid': 'SARS-CoV-2',
    'sars_cov2_rdrp_mn996531': 'SARS-CoV-2',
    'sars_cov1_spike_ay291451': 'SARS-CoV-1',
    'mers_spike_kc881005': 'MERS-CoV',
    'hcov_229e_spike': 'Human coronavirus 229E',
    'hcov_oc43_spike': 'Human coronavirus OC43',
    'hcov_hku1_spike': 'Human coronavirus HKU1',
    'hcov_nl63_spike': 'Human coronavirus NL63',
    
    # Filoviruses
    'ebola_vp35_km034562': 'Ebolavirus',
    'ebola_vp35_kc242800': 'Ebolavirus',
    'ebola_sudan_ay729654': 'Ebolavirus',
    'ebola_bundibugyo': 'Ebolavirus',
    'bundibugyo': 'Ebolavirus',
    'marburg_vp40': 'Marburgvirus',
    'marburg_gp': 'Marburgvirus',
    'reston_vp40': 'Reston ebolavirus',
    'tai_forest': 'Taï Forest virus',
    
    # Arenaviruses
    'lassa_np_ay859705': 'Lassa virus',
    'lassa_gp_ay859704': 'Lassa virus',
    'lassa_l_jq860605': 'Lassa virus',
    'junin_np': 'Junín virus',
    'machupo_np': 'Machupo virus',
    'guanarito_np': 'Guanarito virus',
    'sabia_np': 'Sabia virus',
    'whitewater_np': 'Whitewater Arroyo virus',
    'lcmv_np': 'Lymphocytic choriomeningitis virus',
    
    # Paramyxoviruses
    'measles_f_ab470224': 'Measles virus',
    'measles_m': 'Measles virus',
    'measles_l': 'Measles virus',
    'mumps_f': 'Mumps virus',
    'mumps_np': 'Mumps virus',
    'rsv_g': 'Respiratory syncytial virus',
    'parainfluenza_f': 'Parainfluenza virus',
    'sendai_f': 'Sendai virus',
}

print("="*80)
print("CODON USAGE TABLE DOWNLOAD GUIDE")
print("="*80)
print("\nSince Kazusa doesn't have an easy API, you'll need to download manually.")
print("Here's the organized list:\n")

# Organize by organism
organisms_set = set(organisms_to_download.values())
print(f"Found {len(organisms_set)} unique organisms to download:\n")

for organism in sorted(organisms_set):
    filenames = [k for k, v in organisms_to_download.items() if v == organism]
    print(f"\n{organism}:")
    print(f"  Use for: {', '.join(filenames[:3])}" + ("..." if len(filenames) > 3 else ""))

print("\n" + "="*80)
print("MANUAL DOWNLOAD STEPS:")
print("="*80)
print("""
1. Go to: https://www.ncbi.nlm.nih.gov/codon/table/
2. Search for each unique organism above
3. Download the codon usage table as .txt
4. Save to: data/codon_usage/[organism_name].txt

Example:
- Search: "Influenza A virus"
- Download table
- Save as: data/codon_usage/influenza_a_virus.txt

Then update extract_codon_usage.py to map all filenames to organism codon tables.
""")

print("\nUnique organisms to download codon tables for:")
for i, organism in enumerate(sorted(organisms_set), 1):
    print(f"{i:2d}. {organism}")