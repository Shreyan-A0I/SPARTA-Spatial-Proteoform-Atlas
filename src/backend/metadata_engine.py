"""
Metadata Engine (Phase A: The Blueprint)
=========================================
Constructs a theoretical mass library with structural cross-references.

Key Functions:
- Mass calculation with methionine loss logic
- UniProt/PDB hierarchy management
- Structural data integration (PDB > AlphaFold priority)
"""

import numpy as np
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass


@dataclass
class ProteinMetadata:
    """Stores protein metadata including masses and structural references."""
    uniprot_id: str
    sequence: str
    theoretical_mass: float  # Full mass including N-terminal Met
    theoretical_mass_no_met: float  # Mass after N-terminal Met cleavage
    pdb_ids: List[str]  # PDB structures (X-ray/Cryo-EM)
    alphafold_id: Optional[str]  # AlphaFold prediction ID
    hierarchy_priority: str  # "PDB" or "AlphaFold"


class MetadataEngine:
    """
    The Metadata Engine: Constructs theoretical mass library with structural cross-references.
    
    Mathematical Foundations:
    - Theoretical mass: m_th = Σ(Residue Masses) + 18.01056 Da (H2O)
    - Methionine loss: m_th - 131.0405 Da (if N-terminal Met is cleaved)
    - Hierarchy: PDB (X-ray/Cryo-EM) > AlphaFold (Predicted)
    """
    
    # Amino acid masses (monoisotopic)
    AA_MASSES = {
        'A': 71.03711, 'R': 156.10111, 'N': 114.04293, 'D': 115.02694,
        'C': 103.00919, 'E': 129.04259, 'Q': 128.05858, 'G': 57.02146,
        'H': 137.05891, 'I': 113.08406, 'L': 113.08406, 'K': 128.09496,
        'M': 131.04049, 'F': 147.06841, 'P': 97.05276, 'S': 87.03203,
        'T': 101.04768, 'W': 186.07931, 'Y': 163.06333, 'V': 99.06841
    }
    
    WATER_MASS = 18.01056  # H2O mass
    METHIONINE_MASS = 131.04049  # N-terminal Met mass
    
    def __init__(self):
        """Initialize the metadata engine."""
        self.protein_library: Dict[str, ProteinMetadata] = {}
    
    def calculate_theoretical_mass(self, sequence: str) -> Tuple[float, float]:
        """
        Calculate theoretical mass of a protein sequence.
        
        Formula: m_th = Σ(Residue Masses) + 18.01056 Da
        
        Args:
            sequence: Protein amino acid sequence (single letter code)
            
        Returns:
            Tuple of (full_mass, mass_without_n_term_met)
        """
        # Calculate sum of residue masses
        total_mass = sum(self.AA_MASSES.get(aa, 0) for aa in sequence.upper())
        
        # Add water molecule
        full_mass = total_mass + self.WATER_MASS
        
        # Calculate mass without N-terminal Met (if first residue is Met)
        mass_no_met = full_mass
        if sequence and sequence[0].upper() == 'M':
            mass_no_met = full_mass - self.METHIONINE_MASS
        
        return full_mass, mass_no_met
    
    def fetch_uniprot_data(self, uniprot_id: str) -> Optional[Dict]:
        """
        Fetch protein data from UniProt database.
        
        TODO (Task T1): Implement UniProt API integration
        - Query UniProt REST API
        - Extract sequence, PDB cross-references
        - Handle hierarchy logic (PDB > AlphaFold)
        
        Args:
            uniprot_id: UniProt accession number
            
        Returns:
            Dictionary with protein data or None if not found
        """
        # Placeholder implementation
        # TODO: Implement actual UniProt API calls
        return None
    
    def fetch_pdb_references(self, uniprot_id: str) -> List[str]:
        """
        Fetch PDB structure IDs associated with a UniProt entry.
        
        TODO (Task T1): Implement PDB cross-reference extraction
        - Parse UniProt XML/JSON response
        - Extract PDB IDs from cross-references
        - Filter by experimental method (X-ray, Cryo-EM)
        
        Args:
            uniprot_id: UniProt accession number
            
        Returns:
            List of PDB IDs
        """
        # Placeholder implementation
        return []
    
    def fetch_alphafold_id(self, uniprot_id: str) -> Optional[str]:
        """
        Fetch AlphaFold structure ID for a UniProt entry.
        
        TODO (Task T1): Implement AlphaFold DB lookup
        - Query AlphaFold database
        - Return structure ID if available
        
        Args:
            uniprot_id: UniProt accession number
            
        Returns:
            AlphaFold structure ID or None
        """
        # Placeholder implementation
        return None
    
    def add_protein(self, uniprot_id: str, sequence: str, 
                   pdb_ids: List[str] = None, alphafold_id: str = None) -> ProteinMetadata:
        """
        Add a protein to the metadata library.
        
        Args:
            uniprot_id: UniProt accession number
            sequence: Protein amino acid sequence
            pdb_ids: List of PDB structure IDs
            alphafold_id: AlphaFold structure ID
            
        Returns:
            ProteinMetadata object
        """
        full_mass, mass_no_met = self.calculate_theoretical_mass(sequence)
        
        # Determine hierarchy priority
        hierarchy = "PDB" if pdb_ids else "AlphaFold" if alphafold_id else "None"
        
        metadata = ProteinMetadata(
            uniprot_id=uniprot_id,
            sequence=sequence,
            theoretical_mass=full_mass,
            theoretical_mass_no_met=mass_no_met,
            pdb_ids=pdb_ids or [],
            alphafold_id=alphafold_id,
            hierarchy_priority=hierarchy
        )
        
        self.protein_library[uniprot_id] = metadata
        return metadata
    
    def get_protein_metadata(self, uniprot_id: str) -> Optional[ProteinMetadata]:
        """
        Retrieve protein metadata from the library.
        
        Args:
            uniprot_id: UniProt accession number
            
        Returns:
            ProteinMetadata object or None if not found
        """
        return self.protein_library.get(uniprot_id)
    
    def search_by_mass(self, observed_mass: float, ppm_tolerance: float = 50.0) -> List[ProteinMetadata]:
        """
        Search for proteins matching an observed mass within ppm tolerance.
        
        Args:
            observed_mass: Observed mass from MSI data
            ppm_tolerance: Mass accuracy tolerance in parts per million
            
        Returns:
            List of matching ProteinMetadata objects
        """
        matches = []
        for metadata in self.protein_library.values():
            # Check both masses (with and without N-terminal Met)
            for mass in [metadata.theoretical_mass, metadata.theoretical_mass_no_met]:
                ppm_error = abs(observed_mass - mass) / mass * 1e6
                if ppm_error <= ppm_tolerance:
                    matches.append(metadata)
                    break
        return matches

