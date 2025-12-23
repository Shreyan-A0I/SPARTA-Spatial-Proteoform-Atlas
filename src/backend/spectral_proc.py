"""
Spectral Pipeline (Phase B: The Evidence)
===========================================
Denoises raw MSI data and normalizes across the tissue slice.

Key Functions:
- TIC Normalization
- Peak Picking (Centroiding) with S/N threshold
- Unsupervised Clustering (K-Means/UMAP)
"""

import numpy as np
from typing import Tuple, List, Optional
from dataclasses import dataclass

# Optional import for imzML support
try:
    import pyimzML.ImzMLReader as imzml
except ImportError:
    imzml = None


@dataclass
class SpectralPeak:
    """Represents a detected peak in the mass spectrum."""
    mz: float  # Mass-to-charge ratio
    intensity: float  # Peak intensity
    snr: float  # Signal-to-noise ratio
    pixel_coords: Tuple[int, int]  # Spatial coordinates


@dataclass
class ProcessedSpectrum:
    """Processed spectrum data for a single pixel."""
    pixel_coords: Tuple[int, int]
    mz_array: np.ndarray
    intensity_array: np.ndarray
    tic: float  # Total Ion Current
    peaks: List[SpectralPeak]


class SpectralProcessor:
    """
    The Spectral Pipeline: Denoises raw MSI data and normalizes across tissue slice.
    
    Mathematical Foundations:
    - TIC Normalization: I_norm = I_i / Σ(I_j)
    - Peak Picking: S/N = (I_peak - μ_noise) / σ_noise ≥ 3
    - Cosine Similarity: similarity = (A·B) / (||A|| ||B||)
    """
    
    def __init__(self, snr_threshold: float = 3.0):
        """
        Initialize the spectral processor.
        
        Args:
            snr_threshold: Minimum signal-to-noise ratio for peak detection
        """
        self.snr_threshold = snr_threshold
    
    def load_imzml(self, imzml_path: str, ibd_path: str = None):
        """
        Load imzML file and associated .ibd file.
        
        TODO (Task T2): Implement full imzML loading
        - Handle both continuous and processed imzML formats
        - Load binary data from .ibd file
        - Parse metadata (pixel dimensions, m/z range, etc.)
        
        Args:
            imzml_path: Path to .imzML file
            ibd_path: Path to .ibd file (optional, auto-detected if None)
            
        Returns:
            ImzMLReader object
        """
        if imzml is None:
            raise ImportError("pyimzML is not installed. Install it with: pip install pyimzML")
        # Placeholder implementation
        # TODO: Implement actual imzML loading
        # reader = imzml.ImzMLReader(imzml_path, ibd_path)
        # return reader
        raise NotImplementedError("imzML loading not yet implemented (Task T2)")
    
    def calculate_tic(self, intensity_array: np.ndarray) -> float:
        """
        Calculate Total Ion Current for a spectrum.
        
        Args:
            intensity_array: Array of intensity values
            
        Returns:
            Total Ion Current value
        """
        return np.sum(intensity_array)
    
    def normalize_tic(self, intensity_array: np.ndarray) -> np.ndarray:
        """
        Normalize intensities by Total Ion Current.
        
        Formula: I_norm = I_i / Σ(I_j)
        
        Args:
            intensity_array: Raw intensity array
            
        Returns:
            TIC-normalized intensity array
        """
        tic = self.calculate_tic(intensity_array)
        if tic == 0:
            return intensity_array
        return intensity_array / tic
    
    def estimate_noise(self, intensity_array: np.ndarray, 
                      percentile: float = 10.0) -> Tuple[float, float]:
        """
        Estimate noise statistics from intensity distribution.
        
        Uses lower percentile of intensities as noise baseline.
        
        Args:
            intensity_array: Array of intensity values
            percentile: Percentile to use for noise estimation
            
        Returns:
            Tuple of (mean_noise, std_noise)
        """
        noise_threshold = np.percentile(intensity_array, percentile)
        noise_values = intensity_array[intensity_array <= noise_threshold]
        
        if len(noise_values) == 0:
            return 0.0, 0.0
        
        mean_noise = np.mean(noise_values)
        std_noise = np.std(noise_values)
        
        return mean_noise, std_noise
    
    def calculate_snr(self, peak_intensity: float, mean_noise: float, 
                     std_noise: float) -> float:
        """
        Calculate Signal-to-Noise Ratio for a peak.
        
        Formula: S/N = (I_peak - μ_noise) / σ_noise
        
        Args:
            peak_intensity: Intensity of the peak
            mean_noise: Mean noise level
            std_noise: Standard deviation of noise
            
        Returns:
            Signal-to-noise ratio
        """
        if std_noise == 0:
            return float('inf') if peak_intensity > mean_noise else 0.0
        return (peak_intensity - mean_noise) / std_noise
    
    def peak_picking(self, mz_array: np.ndarray, intensity_array: np.ndarray,
                    pixel_coords: Tuple[int, int]) -> List[SpectralPeak]:
        """
        Detect peaks in a spectrum using S/N threshold.
        
        A peak is kept if: S/N = (I_peak - μ_noise) / σ_noise ≥ threshold
        
        TODO (Task T2): Implement advanced peak picking
        - Local maxima detection
        - Peak width estimation
        - Isotope pattern recognition
        
        Args:
            mz_array: Mass-to-charge ratio array
            intensity_array: Intensity array (should be normalized)
            pixel_coords: Spatial coordinates of the pixel
            
        Returns:
            List of detected SpectralPeak objects
        """
        if len(mz_array) == 0 or len(intensity_array) == 0:
            return []
        
        # Estimate noise
        mean_noise, std_noise = self.estimate_noise(intensity_array)
        
        # Find local maxima (simple implementation)
        peaks = []
        for i in range(1, len(intensity_array) - 1):
            if (intensity_array[i] > intensity_array[i-1] and 
                intensity_array[i] > intensity_array[i+1]):
                
                # Calculate S/N
                snr = self.calculate_snr(intensity_array[i], mean_noise, std_noise)
                
                # Keep peak if S/N exceeds threshold
                if snr >= self.snr_threshold:
                    peaks.append(SpectralPeak(
                        mz=mz_array[i],
                        intensity=intensity_array[i],
                        snr=snr,
                        pixel_coords=pixel_coords
                    ))
        
        return peaks
    
    def process_spectrum(self, mz_array: np.ndarray, intensity_array: np.ndarray,
                        pixel_coords: Tuple[int, int]) -> ProcessedSpectrum:
        """
        Process a single spectrum: normalize and detect peaks.
        
        Args:
            mz_array: Mass-to-charge ratio array
            intensity_array: Raw intensity array
            pixel_coords: Spatial coordinates
            
        Returns:
            ProcessedSpectrum object
        """
        # Normalize by TIC
        normalized_intensity = self.normalize_tic(intensity_array)
        
        # Calculate TIC
        tic = self.calculate_tic(intensity_array)
        
        # Detect peaks
        peaks = self.peak_picking(mz_array, normalized_intensity, pixel_coords)
        
        return ProcessedSpectrum(
            pixel_coords=pixel_coords,
            mz_array=mz_array,
            intensity_array=normalized_intensity,
            tic=tic,
            peaks=peaks
        )
    
    def cosine_similarity(self, spectrum1: np.ndarray, spectrum2: np.ndarray) -> float:
        """
        Calculate cosine similarity between two spectra.
        
        Formula: similarity = (A·B) / (||A|| ||B||)
        
        Args:
            spectrum1: First intensity array
            spectrum2: Second intensity array (must be same length)
            
        Returns:
            Cosine similarity value (0-1)
        """
        if len(spectrum1) != len(spectrum2):
            raise ValueError("Spectra must have the same length")
        
        dot_product = np.dot(spectrum1, spectrum2)
        norm1 = np.linalg.norm(spectrum1)
        norm2 = np.linalg.norm(spectrum2)
        
        if norm1 == 0 or norm2 == 0:
            return 0.0
        
        return dot_product / (norm1 * norm2)
    
    def cluster_pixels(self, processed_spectra: List[ProcessedSpectrum],
                      method: str = "kmeans", n_clusters: int = 5) -> np.ndarray:
        """
        Cluster pixels based on spectral similarity.
        
        TODO (Task T2): Implement clustering algorithms
        - K-Means clustering
        - UMAP dimensionality reduction + clustering
        - Use cosine similarity for distance metric
        
        Args:
            processed_spectra: List of processed spectra
            method: Clustering method ("kmeans" or "umap")
            n_clusters: Number of clusters (for K-Means)
            
        Returns:
            Array of cluster labels for each pixel
        """
        # Placeholder implementation
        # TODO: Implement actual clustering
        return np.zeros(len(processed_spectra), dtype=int)

