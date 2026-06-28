"""
Download viral sequences from NCBI GenBank in batch
Requires: biopython, requests
"""

from Bio import Entrez
import os

# Set your email (NCBI requires this)
Entrez.email = "didenkoolena2504@gmail.com"

# List of sequences to download
sequences_to_download = {
    # Influenza
    'influenza_a_ha_cy163680': 'CY163680.1',
    'influenza_a_ha_jq708657': 'JQ708657.1',
    'influenza_a_ha_kf821935': 'KF821935.1',
    'influenza_a_na_cy163681': 'CY163681.1',
    'influenza_a_pa_cy163684': 'CY163684.1',
    'influenza_a_m1_cy163686': 'CY163686.1',
    'influenza_b_ha_lc107905': 'LC107905.1',
    'influenza_c_ha_ab211639': 'AB211639.1',
    'h1n1_ha_gq152857': 'GQ152857.1',
    
    # Coronaviruses
    'sars_cov2_spike_ov109917': 'OV109917.1',
    'sars_cov2_nucleocapsid': 'NC_045512.2',
    'sars_cov2_rdrp_mn996531': 'MN996531.1',
    'sars_cov1_spike_ay291451': 'AY291451.1',
    'mers_spike_kc881005': 'KC881005.1',
    'hcov_229e_spike': 'AY585882.1',
    'hcov_oc43_spike': 'AY585887.1',
    'hcov_hku1_spike': 'AY500324.1',
    'hcov_nl63_spike': 'AY567487.1',
    
    # Filoviruses
    'ebola_vp35_km034562': 'KM034562.1',
    'ebola_vp35_kc242800': 'KC242800.1',
    'ebola_sudan_ay729654': 'AY729654.1',
    'ebola_bundibugyo': 'FJ217161.1',
    'marburg_vp40': 'NC_001608.1',
    'marburg_gp': 'NC_001608.1',
    'reston_vp40': 'NC_004161.1',
    'tai_forest': 'NC_014372.1',
    'bundibugyo': 'NC_014373.1',
    
    # Arenaviruses
    'lassa_np_ay859705': 'AY859705.1',
    'lassa_gp_ay859704': 'AY859704.1',
    'lassa_l_jq860605': 'JQ860605.1',
    'junin_np': 'NC_005081.1',
    'machupo_np': 'NC_005077.1',
    'guanarito_np': 'NC_005078.1',
    'sabia_np': 'NC_005080.1',
    'whitewater_np': 'NC_006437.1',
    'lcmv_np': 'NC_001423.1',
    
    # Paramyxoviruses
    'measles_f_ab470224': 'AB470224.1',
    'measles_m': 'M81897.1',
    'measles_l': 'AB470225.1',
    'mumps_f': 'AB010145.1',
    'mumps_np': 'M33221.1',
    'rsv_f': 'M17216.1',
    'rsv_g': 'M17211.1',
    'parainfluenza_f': 'M29391.1',
    'sendai_f': 'M58450.1',
}

def download_sequence(accession, filename):
    """Download sequence from NCBI and save to FASTA"""
    try:
        print(f"Downloading {accession}...", end=' ')
        handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
        record = handle.read()
        
        filepath = os.path.join('data/reference_sequences', f"{filename}.fasta")
        with open(filepath, 'w') as f:
            f.write(record)
        
        print(f"✅ Saved to {filepath}")
        return True
    
    except Exception as e:
        print(f"❌ Error: {e}")
        return False

# Main
if __name__ == '__main__':
    print("="*80)
    print("BATCH DOWNLOADING VIRAL SEQUENCES FROM NCBI")
    print("="*80 + "\n")
    
    os.makedirs('data/reference_sequences', exist_ok=True)
    
    downloaded = 0
    failed = 0
    
    for name, accession in sequences_to_download.items():
        if download_sequence(accession, name):
            downloaded += 1
        else:
            failed += 1
    
    print("\n" + "="*80)
    print(f"Downloaded: {downloaded}/{len(sequences_to_download)}")
    print(f"Failed: {failed}/{len(sequences_to_download)}")
    print("="*80 + "\n")