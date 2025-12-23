"""
Backend modules for SPARTA.
"""

from .metadata_engine import MetadataEngine, ProteinMetadata
from .spectral_proc import SpectralProcessor, ProcessedSpectrum, SpectralPeak
from .discovery_logic import DiscoveryMatcher, ProteoformMatch, PTMDefinition

__all__ = [
    'MetadataEngine',
    'ProteinMetadata',
    'SpectralProcessor',
    'ProcessedSpectrum',
    'SpectralPeak',
    'DiscoveryMatcher',
    'ProteoformMatch',
    'PTMDefinition',
]

