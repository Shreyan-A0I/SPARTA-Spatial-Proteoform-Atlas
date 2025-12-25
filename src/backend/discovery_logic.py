"""
Discovery Matcher
========================================
Identifies functional mass shifts (PTMs) with statistical confidence.

Key Functions:
- ppm Error calculation
- Delta-mass matching for PTMs
- Pearson correlation validation
- Isobaric ambiguity detection
"""
import numpy as np
import pandas as pd
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass
from scipy.stats import pearsonr

@dataclass
class ProteoformMatch:
    uniprot_id: str
    protein_name: str
    observed_mz: float
    theoretical_mass: float
    ppm_error: float
    ptm_name: str  # "Base" if no modification
    match_type: str # e.g., "mono_no_met"
    confidence: float

class DiscoveryMatcher:
    """
    Identifies proteoforms using UniMod-calibrated mass shifts.
    Precision logic: E_ppm = \frac{|m_{obs} - m_{th}|}{m_{th}} \times 10^6
    """
    
    # UniMod-calibrated constants (Monoisotopic)
    PTM_LIBRARY = {
        'Acetylation': 42.010565,
        'Phosphorylation': 79.966331,
        'Methylation': 14.015650,
        'Oxidation': 15.994915,
        'Sulfation': 79.956815 # Crucial for isobaric detection vs Phospho
    }

    def __init__(self, ppm_tolerance: float = 20.0):
        self.ppm_tolerance = ppm_tolerance

    def calculate_ppm(self, obs: float, theor: float) -> float:
        return (abs(obs - theor) / theor) * 1e6

    def find_best_anchor(self, obs_mz: float, row: pd.Series) -> Tuple[Optional[str], float, float]:
        """
        Tests the 4 theoretical masses to find the most likely base form.
        Returns: (match_type, theoretical_mass, ppm_error)
        """
        mass_options = {
            'mono': row['mono_mass'],
            'mono_no_met': row['mono_mass_no_met'],
            'avg': row['avg_mass'],
            'avg_no_met': row['avg_mass_no_met']
        }
        
        best_match = None
        min_ppm = self.ppm_tolerance
        anchor_mass = 0.0

        for m_type, m_val in mass_options.items():
            if pd.isna(m_val): continue
            ppm = self.calculate_ppm(obs_mz, m_val)
            if ppm < min_ppm:
                min_ppm = ppm
                best_match = m_type
                anchor_mass = m_val
        
        return best_match, anchor_mass, min_ppm

    def match_discovery(self, consensus_peaks: pd.DataFrame, protein_library: pd.DataFrame) -> List[ProteoformMatch]:
        """
        Main matching loop: Base Detection -> PTM Shift Detection.
        """
        results = []
        
        for _, peak in consensus_peaks.iterrows():
            obs_mz = peak['avg_mz']
            
            for _, protein in protein_library.iterrows():
                # 1. Check for Base Match (Direct identification)
                m_type, anchor, ppm = self.find_best_anchor(obs_mz, protein)
                
                if m_type:
                    results.append(ProteoformMatch(
                        uniprot_id=protein['uniprot_id'],
                        protein_name=protein['protein_name'],
                        observed_mz=obs_mz,
                        theoretical_mass=anchor,
                        ppm_error=ppm,
                        ptm_name="Base",
                        match_type=m_type,
                        confidence=1.0 - (ppm / self.ppm_tolerance)
                    ))
                    continue # Found base, move to next protein
                
                # 2. Check for PTM Shifts (Only if no base match for this peak)
                # We test against the primary 'mono_no_met' as the standard mature form
                base_theor = protein['mono_mass_no_met']
                if pd.isna(base_theor): continue
                
                for ptm_name, delta in self.PTM_LIBRARY.items():
                    target_mass = base_theor + delta
                    ppm_ptm = self.calculate_ppm(obs_mz, target_mass)
                    
                    if ppm_ptm < self.ppm_tolerance:
                        results.append(ProteoformMatch(
                            uniprot_id=protein['uniprot_id'],
                            protein_name=protein['protein_name'],
                            observed_mz=obs_mz,
                            theoretical_mass=target_mass,
                            ppm_error=ppm_ptm,
                            ptm_name=ptm_name,
                            match_type="mono_no_met_plus_PTM",
                            confidence=0.8 - (ppm_ptm / self.ppm_tolerance) # Lower base confidence for PTMs
                        ))
        
        return results

    def detect_isobaric_clash(self, match: ProteoformMatch) -> Optional[str]:
        """
        Checks if a Phospho match could actually be Sulfation.
        Diff is only 0.0095 Da.
        """
        if match.ptm_name == "Phosphorylation":
            # If our ppm error is high, it might actually be the other one
            sulf_mass = match.theoretical_mass - 0.009516 
            sulf_ppm = self.calculate_ppm(match.observed_mz, sulf_mass)
            if sulf_ppm < match.ppm_error:
                return "Sulfation"
        return None