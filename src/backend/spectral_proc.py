"""
Spectral Pipeline
===========================================
Denoises raw MSI data and normalizes across the tissue slice.

Key Functions:
- TIC Normalization
- Peak Picking (Centroiding) with S/N threshold
- Mass Consensus: Aggregates peaks across pixels to identify significant features
- Excel Export: Generates a detailed report of the spectral analysis
"""

import numpy as np
import pandas as pd
from typing import Tuple, List, Dict
from dataclasses import dataclass, asdict
from pyimzml.ImzMLParser import ImzMLParser
from scipy.signal import find_peaks

@dataclass
class SpectralPeak:
    """Represents a detected peak in a single pixel."""
    mz: float
    intensity: float
    snr: float
    pixel_coords: Tuple[int, int]

class SpectralProcessor:
    """
    SPARTA Spectral Pipeline: Pre-processing and Evidence Extraction.
    Handles TIC normalization, peak picking, and mass consensus.
    """
    
    def __init__(self, snr_threshold: float = 3.0):
        self.snr_threshold = snr_threshold
        self.parser = None

    def load_data(self, imzml_path: str):
        """Loads MSI data into the parser."""
        self.parser = ImzMLParser(imzml_path)
        print(f"✅ Loaded {len(self.parser.coordinates)} pixels from {imzml_path}")

    def get_tic_map(self) -> np.ndarray:
        """Generates a 2D intensity heatmap (Tissue Silhouette)."""
        coords = np.array(self.parser.coordinates)
        xs, ys = coords[:, 0], coords[:, 1]
        x_max, y_max = xs.max(), ys.max()
        tic_map = np.zeros((y_max + 1, x_max + 1), dtype=np.float32)

        for i, (x, y, *_) in enumerate(self.parser.coordinates):
            _, intensities = self.parser.getspectrum(i)
            tic_map[y, x] = np.sum(intensities)
        return tic_map

    def process_and_extract(self) -> List[SpectralPeak]:
        if self.parser is None:
            raise RuntimeError("Call load_data() before process_and_extract().")

        master_peaks: List[SpectralPeak] = []
        print("🚀 Processing pixels and picking peaks...")

        # Diagnostic: Let's look at the first pixel to see what we're dealing with
        mzs, intensities = self.parser.getspectrum(0)
        print(f"DEBUG: First pixel max raw intensity: {np.max(intensities)}")

        for i, coords in enumerate(self.parser.coordinates):
            mzs, intensities = self.parser.getspectrum(i)

            tic = np.sum(intensities)
            if tic <= 0: continue
            
            norm_intensities = intensities / tic

            # Use a slightly higher percentile if 10th is too low (e.g., 25th)
            noise_floor = np.percentile(norm_intensities, 25)
            noise_mask = norm_intensities <= noise_floor
            
            # Fallback if noise estimation is impossible
            std_noise = np.std(norm_intensities[noise_mask])
            if std_noise <= 0:
                # If spectrum is too sparse, use a global small value as noise baseline
                std_noise = np.mean(norm_intensities) * 0.1 

            # Peak Picking
            # If find_peaks is too strict, we lower the bar
            snr_min_height = noise_floor + (self.snr_threshold * std_noise)
            indices, _ = find_peaks(norm_intensities, height=snr_min_height)

            for idx in indices:
                p_intensity = norm_intensities[idx]
                snr = (p_intensity - noise_floor) / (std_noise + 1e-9)

                master_peaks.append(SpectralPeak(
                    mz=float(mzs[idx]),
                    intensity=float(p_intensity),
                    snr=float(snr),
                    pixel_coords=(coords[0], coords[1])
                ))
            
            # Print progress every 5000 pixels
            if i % 5000 == 0 and i > 0:
                print(f"   ...processed {i} pixels ({len(master_peaks)} peaks found so far)")
        
        print(f"🎯 Extracted {len(master_peaks)} raw peak detections.")
        return master_peaks

    def export_consensus_to_excel(self, master_peaks: List[SpectralPeak], 
                                 filename: str = "SPARTA_Spectral_Sanity_Check.xlsx", 
                                 tolerance_da: float = 0.5):
        # PROTECT AGAINST EMPTY LIST
        if not master_peaks:
            print("⚠️ No peaks were detected. Skipping Excel export. Try lowering snr_threshold.")
            return

        print(f"📊 Binning peaks...")
        df = pd.DataFrame([asdict(p) for p in master_peaks])
        
        # Ensure 'mz' column exists (it should now if master_peaks is not empty)
        if 'mz' not in df.columns:
            print("⚠️ Data error: 'mz' column missing in peak list.")
            return

        df['mz_bin'] = np.round(df['mz'] / tolerance_da) * tolerance_da
        
        # Aggregate to find Consensus Peaks
        total_pixels = len(self.parser.coordinates)
        summary = df.groupby('mz_bin').agg({
            'mz': 'mean',           # Actual average mz in that bin
            'intensity': 'mean',    # Mean normalized intensity
            'snr': 'max',           # Maximum SNR observed
            'pixel_coords': 'count' # Spatial frequency
        }).rename(columns={'pixel_coords': 'spatial_count', 'mz': 'avg_mz'})

        # Filter: Must be in >1% of tissue to avoid exporting random noise
        min_pixels = max(5, int(total_pixels * 0.01))
        consensus_df = summary[summary['spatial_count'] >= min_pixels].copy()
        consensus_df = consensus_df.sort_values(by='spatial_count', ascending=False)

        # Export
        with pd.ExcelWriter(filename) as writer:
            consensus_df.to_excel(writer, sheet_name='Consensus Peaks')
            
            # Metadata Sheet
            metadata = pd.DataFrame({
                "Metric": ["Total Pixels", "Total Raw Peaks", "Consensus Peaks (>1%)", "SNR Threshold"],
                "Value": [total_pixels, len(master_peaks), len(consensus_df), self.snr_threshold]
            })
            metadata.to_excel(writer, sheet_name='Run Summary', index=False)

        print(f"✨ Excel report generated: {filename}")

# --- EXAMPLE EXECUTION ---
if __name__ == "__main__":
    processor = SpectralProcessor(snr_threshold=3.5)
    
    # Update these paths to your actual local files
    try:
        processor.load_data("data/mouse.imzML")
        raw_peaks = processor.process_and_extract()
        processor.export_consensus_to_excel(raw_peaks)
    except FileNotFoundError:
        print("❌ File not found. Please ensure mouse.imzML is in the data/ folder.")