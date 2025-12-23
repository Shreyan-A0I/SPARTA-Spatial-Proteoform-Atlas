"""
Example usage of SPARTA modules.

This script demonstrates how to use the core components of SPARTA
for spatial proteoform analysis.
"""

import sys
import os

# Add src to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

from backend.metadata_engine import MetadataEngine
from backend.spectral_proc import SpectralProcessor
from backend.discovery_logic import DiscoveryMatcher
from utils.sasa_engine import SASAEngine
import numpy as np


def example_metadata_engine():
    """Example: Using the Metadata Engine."""
    print("=" * 60)
    print("Example 1: Metadata Engine")
    print("=" * 60)
    
    engine = MetadataEngine()
    
    # Example protein sequence (insulin)
    sequence = "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN"
    
    # Calculate theoretical mass
    full_mass, mass_no_met = engine.calculate_theoretical_mass(sequence)
    
    print(f"Sequence length: {len(sequence)}")
    print(f"Theoretical mass (with N-term Met): {full_mass:.2f} Da")
    print(f"Theoretical mass (without N-term Met): {mass_no_met:.2f} Da")
    print(f"Mass difference: {full_mass - mass_no_met:.2f} Da")
    print()


def example_spectral_processing():
    """Example: Using the Spectral Processor."""
    print("=" * 60)
    print("Example 2: Spectral Processing")
    print("=" * 60)
    
    processor = SpectralProcessor(snr_threshold=3.0)
    
    # Generate sample spectrum data
    mz_array = np.linspace(500, 5000, 1000)
    intensity_array = np.random.rand(1000) * 1000
    # Add some peaks
    intensity_array[200] = 5000
    intensity_array[500] = 8000
    intensity_array[800] = 3000
    
    # Calculate TIC
    tic = processor.calculate_tic(intensity_array)
    print(f"Total Ion Current: {tic:.2f}")
    
    # Normalize by TIC
    normalized = processor.normalize_tic(intensity_array)
    print(f"Normalized intensity range: [{normalized.min():.6f}, {normalized.max():.6f}]")
    
    # Detect peaks
    peaks = processor.peak_picking(mz_array, normalized, (0, 0))
    print(f"Detected {len(peaks)} peaks with S/N >= 3.0")
    for i, peak in enumerate(peaks[:5], 1):  # Show first 5
        print(f"  Peak {i}: m/z={peak.mz:.2f}, intensity={peak.intensity:.6f}, S/N={peak.snr:.2f}")
    print()


def example_discovery_matcher():
    """Example: Using the Discovery Matcher."""
    print("=" * 60)
    print("Example 3: Discovery Matcher")
    print("=" * 60)
    
    matcher = DiscoveryMatcher(ppm_tolerance=50.0)
    
    # Example: Base protein mass
    base_mass = 5000.0  # Da
    
    # Example: Observed mass (base + phosphorylation)
    observed_mass = base_mass + 79.9663  # Phosphorylation
    
    # Calculate ppm error
    ppm_error = matcher.calculate_ppm_error(observed_mass, base_mass + 79.9663)
    print(f"Base mass: {base_mass:.2f} Da")
    print(f"Observed mass: {observed_mass:.2f} Da")
    print(f"PPM error: {ppm_error:.2f} ppm")
    
    # Find PTM candidates
    candidates = matcher.find_ptm_candidates(observed_mass, base_mass)
    print(f"\nFound {len(candidates)} PTM candidates:")
    for ptm_name, delta, ppm in candidates:
        print(f"  {ptm_name}: Δ={delta:.4f} Da, error={ppm:.2f} ppm")
    
    # Check for isobaric ambiguity
    ambiguous = matcher.detect_isobaric_ambiguity(observed_mass, base_mass)
    if ambiguous:
        print(f"\n⚠️  Isobaric ambiguity detected: {', '.join(ambiguous)}")
    print()


def example_sasa_engine():
    """Example: Using the SASA Engine."""
    print("=" * 60)
    print("Example 4: SASA Engine")
    print("=" * 60)
    
    sasa = SASAEngine()
    
    # Check accessibility threshold
    test_sasa_values = [10.0, 25.0, 50.0, 100.0]
    print(f"SASA Threshold: {sasa.SASA_THRESHOLD} Å²")
    print("\nAccessibility check:")
    for sasa_value in test_sasa_values:
        is_accessible = sasa.check_residue_accessibility(sasa_value)
        status = "✓ Accessible" if is_accessible else "✗ Not accessible"
        print(f"  SASA = {sasa_value:6.1f} Å²: {status}")
    print()


def main():
    """Run all examples."""
    print("\n" + "=" * 60)
    print("SPARTA: Spatial Proteoform Atlas - Example Usage")
    print("=" * 60 + "\n")
    
    try:
        example_metadata_engine()
        example_spectral_processing()
        example_discovery_matcher()
        example_sasa_engine()
        
        print("=" * 60)
        print("All examples completed successfully!")
        print("=" * 60)
        print("\nNote: Some features require full implementation (see TODO comments in code)")
        print("      and external data sources (UniProt, PDB, imzML files).")
        
    except Exception as e:
        print(f"\nError running examples: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()

