"""
Discovery Matcher (Phase C: The Bridge)
========================================
Identifies functional mass shifts (PTMs) with statistical confidence.

Key Functions:
- ppm Error calculation
- Delta-mass matching for PTMs
- Pearson correlation validation
- Isobaric ambiguity detection
"""

import numpy as np
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass
from scipy.stats import pearsonr


@dataclass
class PTMDefinition:
    """Definition of a Post-Translational Modification."""
    name: str
    delta_mass: float  # Mass shift in Da
    target_residues: List[str]  # Amino acids that can be modified


@dataclass
class ProteoformMatch:
    """Represents a matched proteoform with PTM."""
    uniprot_id: str
    base_mass: float
    observed_mass: float
    ppm_error: float
    ptm_name: str
    ptm_delta: float
    spatial_distribution: np.ndarray  # Intensity across pixels
    correlation_score: float  # Pearson r with base protein
    confidence_score: float  # Combined confidence score
    isobaric_ambiguity: List[str]  # List of possible PTMs if ambiguous


class DiscoveryMatcher:
    """
    The Discovery Matcher: Identifies functional mass shifts (PTMs) with statistical confidence.
    
    Mathematical Foundations:
    - ppm Error: E_ppm = |m_obs - m_th| / m_th × 10^6 ≤ 50 ppm
    - Delta Matcher: m_obs ≈ (m_th + Δ_PTM)
    - Pearson Correlation: r = Σ(X_i - X̄)(Y_i - Ȳ) / √(Σ(X_i - X̄)² Σ(Y_i - Ȳ)²)
    - Confidence Threshold: r > 0.6 for high-probability proteoform
    """
    
    # Common PTM mass shifts (monoisotopic)
    PTM_LIBRARY = {
        'Acetylation': PTMDefinition('Acetylation', 42.0106, ['K', 'N-term']),
        'Phosphorylation': PTMDefinition('Phosphorylation', 79.9663, ['S', 'T', 'Y']),
        'Methylation': PTMDefinition('Methylation', 14.0157, ['K', 'R']),
        'Sulfation': PTMDefinition('Sulfation', 79.9568, ['Y', 'T']),
        'Hydroxylation': PTMDefinition('Hydroxylation', 15.9949, ['P', 'K']),
        'Ubiquitination': PTMDefinition('Ubiquitination', 114.0429, ['K']),
    }
    
    DEFAULT_PPM_TOLERANCE = 50.0
    MIN_CORRELATION_THRESHOLD = 0.6
    ISOBARIC_THRESHOLD = 0.01  # Da - mass difference threshold for isobaric ambiguity
    
    def __init__(self, ppm_tolerance: float = DEFAULT_PPM_TOLERANCE):
        """
        Initialize the discovery matcher.
        
        Args:
            ppm_tolerance: Mass accuracy tolerance in parts per million
        """
        self.ppm_tolerance = ppm_tolerance
    
    def calculate_ppm_error(self, observed_mass: float, theoretical_mass: float) -> float:
        """
        Calculate mass error in parts per million.
        
        Formula: E_ppm = |m_obs - m_th| / m_th × 10^6
        
        Args:
            observed_mass: Observed mass from MSI data
            theoretical_mass: Theoretical mass from protein sequence
            
        Returns:
            ppm error value
        """
        if theoretical_mass == 0:
            return float('inf')
        return abs(observed_mass - theoretical_mass) / theoretical_mass * 1e6
    
    def is_valid_match(self, observed_mass: float, theoretical_mass: float) -> bool:
        """
        Check if a mass match is within ppm tolerance.
        
        A match is valid if: E_ppm ≤ tolerance
        
        Args:
            observed_mass: Observed mass
            theoretical_mass: Theoretical mass
            
        Returns:
            True if match is within tolerance
        """
        ppm_error = self.calculate_ppm_error(observed_mass, theoretical_mass)
        return ppm_error <= self.ppm_tolerance
    
    def find_ptm_candidates(self, observed_mass: float, 
                           base_mass: float) -> List[Tuple[str, float, float]]:
        """
        Find PTM candidates by delta-mass matching.
        
        Looks for: m_obs ≈ (m_th + Δ_PTM)
        
        Args:
            observed_mass: Observed mass from MSI
            base_mass: Base protein mass
            
        Returns:
            List of tuples: (ptm_name, delta_mass, ppm_error)
        """
        candidates = []
        mass_difference = observed_mass - base_mass
        
        for ptm_name, ptm_def in self.PTM_LIBRARY.items():
            # Check if mass difference matches PTM delta
            ppm_error = self.calculate_ppm_error(
                observed_mass, 
                base_mass + ptm_def.delta_mass
            )
            
            if ppm_error <= self.ppm_tolerance:
                candidates.append((ptm_name, ptm_def.delta_mass, ppm_error))
        
        return candidates
    
    def detect_isobaric_ambiguity(self, observed_mass: float, 
                                  base_mass: float) -> List[str]:
        """
        Detect isobaric PTM ambiguity (e.g., Phosphorylation vs Sulfation).
        
        Isobaric overlap occurs when mass differences are < 0.01 Da.
        
        Args:
            observed_mass: Observed mass
            base_mass: Base protein mass
            
        Returns:
            List of ambiguous PTM names
        """
        ambiguous = []
        mass_difference = observed_mass - base_mass
        
        # Get all PTM candidates
        candidates = self.find_ptm_candidates(observed_mass, base_mass)
        
        if len(candidates) <= 1:
            return []
        
        # Check for close mass differences
        candidate_masses = [c[1] for c in candidates]
        for i, mass1 in enumerate(candidate_masses):
            for j, mass2 in enumerate(candidate_masses[i+1:], i+1):
                if abs(mass1 - mass2) < self.ISOBARIC_THRESHOLD:
                    ambiguous.extend([candidates[i][0], candidates[j][0]])
        
        return list(set(ambiguous))  # Remove duplicates
    
    def calculate_pearson_correlation(self, distribution1: np.ndarray,
                                     distribution2: np.ndarray) -> Tuple[float, float]:
        """
        Calculate Pearson correlation coefficient between two spatial distributions.
        
        Formula: r = Σ(X_i - X̄)(Y_i - Ȳ) / √(Σ(X_i - X̄)² Σ(Y_i - Ȳ)²)
        
        Args:
            distribution1: First spatial distribution (e.g., modified peak)
            distribution2: Second spatial distribution (e.g., base protein peak)
            
        Returns:
            Tuple of (correlation_coefficient, p_value)
        """
        if len(distribution1) != len(distribution2):
            raise ValueError("Distributions must have the same length")
        
        # Remove NaN and infinite values
        mask = np.isfinite(distribution1) & np.isfinite(distribution2)
        if np.sum(mask) < 2:
            return 0.0, 1.0
        
        correlation, p_value = pearsonr(distribution1[mask], distribution2[mask])
        return correlation, p_value
    
    def validate_proteoform(self, modified_distribution: np.ndarray,
                           base_distribution: np.ndarray) -> Tuple[bool, float]:
        """
        Validate a proteoform by checking spatial correlation.
        
        If a modification is real, the spatial distribution of the "Modified Peak" (X)
        must correlate with its "Base Protein Peak" (Y).
        
        Confidence Threshold: r > 0.6 suggests high-probability proteoform.
        
        Args:
            modified_distribution: Spatial distribution of modified peak
            base_distribution: Spatial distribution of base protein peak
            
        Returns:
            Tuple of (is_valid, correlation_score)
        """
        correlation, p_value = self.calculate_pearson_correlation(
            modified_distribution, base_distribution
        )
        
        is_valid = correlation >= self.MIN_CORRELATION_THRESHOLD
        return is_valid, correlation
    
    def calculate_confidence_score(self, ppm_error: float, correlation: float,
                                  sasa: float = 0.0, weights: Tuple[float, float, float] = None) -> float:
        """
        Calculate combined confidence score for a proteoform discovery.
        
        Formula: C = w₁(E_ppm) + w₂(r) + w₃(SASA)
        
        Args:
            ppm_error: Mass accuracy error in ppm
            correlation: Pearson correlation coefficient
            sasa: Solvent Accessible Surface Area (optional)
            weights: Tuple of weights (w1, w2, w3) for each component
            
        Returns:
            Combined confidence score (0-1, higher is better)
        """
        if weights is None:
            weights = (0.3, 0.5, 0.2)  # Default weights
        
        w1, w2, w3 = weights
        
        # Normalize ppm_error (inverse, lower is better)
        # Assuming max ppm_error of 50, normalize to 0-1
        ppm_score = max(0, 1 - (ppm_error / self.ppm_tolerance))
        
        # Correlation is already 0-1
        corr_score = max(0, correlation)
        
        # Normalize SASA (assuming max SASA of 200 Å²)
        sasa_score = min(1.0, sasa / 200.0) if sasa > 0 else 0.5
        
        confidence = w1 * ppm_score + w2 * corr_score + w3 * sasa_score
        return confidence
    
    def match_proteoforms(self, observed_masses: np.ndarray,
                         base_masses: Dict[str, float],
                         spatial_distributions: Dict[float, np.ndarray]) -> List[ProteoformMatch]:
        """
        Match observed masses to proteoforms with PTMs.
        
        TODO (Task T3/T4): Implement full matching pipeline
        - Iterate through observed masses
        - Find base protein matches
        - Calculate PTM deltas
        - Validate with spatial correlation
        - Calculate confidence scores
        
        Args:
            observed_masses: Array of observed m/z values
            base_masses: Dictionary mapping UniProt IDs to theoretical masses
            spatial_distributions: Dictionary mapping m/z to spatial intensity arrays
            
        Returns:
            List of ProteoformMatch objects
        """
        matches = []
        
        # Placeholder implementation
        # TODO: Implement full matching logic
        
        return matches

