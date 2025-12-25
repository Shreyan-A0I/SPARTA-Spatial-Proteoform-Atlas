"""
Metadata Engine
=========================================
Constructs a theoretical mass library with structural cross-references.

Key Functions:
- Mass calculation with methionine loss logic
- UniProt/PDB hierarchy management
- Structural data integration (PDB > AlphaFold priority)
"""

import requests
import numpy as np
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass, asdict
import pandas as pd


@dataclass
class ProteinMetadata:
    """Stores protein metadata including masses and structural references."""
    uniprot_id: str
    entry_name: str
    protein_name: str
    sequence: str
    # Monoisotopic masses (Calculated from sequence)
    mono_mass: float
    mono_mass_no_met: float
    # Average masses (Fetched from UniProt)
    avg_mass: float
    avg_mass_no_met: float
    
    pdb_ids: List[str]
    alphafold_id: Optional[str]
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
        self.protein_library: Dict[str, ProteinMetadata] = {}
    
        """
        Note: N-terminal Methionine Exopeptidase is an enzyme that often cleaves the 
        starting Methionine if the second residue is small (like Alanine, Serine, or Glycine)
        In a tissue slice, you don't know if the protein exists in its "Full" form or its "Met-cleaved" form
        so to search for proteoforms, both masses are calculated and stored.
        By checking both, you are effectively looking for two potential "Base Peaks."
        """
    def calculate_mono_mass(self, sequence: str) -> Tuple[float, float]:
        """Calculates monoisotopic mass and met-cleaved version."""
        total = sum(self.AA_MASSES.get(aa, 0) for aa in sequence.upper()) + self.WATER_MASS
        no_met = total - self.METHIONINE_MASS if sequence.startswith('M') else total
        return total, no_met
    
    def fetch_proteome(self, taxonomy_id: int, size: int = 500):
        """
        Generalized fetcher for any organism.
        Args:
            taxonomy_id: 10090 (Mouse), 9606 (Human), etc.
            size: Number of reviewed proteins to fetch.
        """
        print(f"Fetching proteome for Taxonomy ID: {taxonomy_id}...")
        
        base_url = "https://rest.uniprot.org/uniprotkb/search"
        # query: Reviewed only, specific taxonomy
        # fields: mass (Avg), sequence, pdb, alphafold
        params = {
            "query": f"taxonomy_id:{taxonomy_id} AND reviewed:true",
            "fields": "accession,id,protein_name,mass,sequence,xref_pdb,xref_alphafolddb",
            "format": "json",
            "size": size,
        }

        try:
            response = requests.get(base_url, params=params)
            response.raise_for_status()
            data = response.json()
            
            for entry in data.get('results', []):
                uid = entry['primaryAccession']
                entry_name = entry['uniProtkbId']
                name = entry['proteinDescription']['recommendedName']['fullName']['value']
                seq = entry['sequence']['value']
                
                # 1. Get Average Mass from UniProt
                avg_f = float(entry['sequence']['molWeight'])
                avg_no_m = avg_f - self.METHIONINE_MASS if seq.startswith('M') else avg_f
                
                # 2. Calculate Monoisotopic Mass
                mono_f, mono_no_m = self.calculate_mono_mass(seq)
                
                # 3. Extract Structure IDs
                xrefs = entry.get('uniProtKBCrossReferences', [])
                pdbs = [x['id'] for x in xrefs if x['database'] == 'PDB']
                af_id = next((x['id'] for x in xrefs if x['database'] == 'AlphaFoldDB'), None)
                
                hierarchy = "PDB" if pdbs else "AlphaFold" if af_id else "None"

                metadata = ProteinMetadata(
                    uniprot_id=uid,
                    entry_name=entry_name,
                    protein_name=name,
                    sequence=seq,
                    mono_mass=mono_f,
                    mono_mass_no_met=mono_no_m,
                    avg_mass=avg_f,
                    avg_mass_no_met=avg_no_m,
                    pdb_ids=pdbs,
                    alphafold_id=af_id,
                    hierarchy_priority=hierarchy
                )
                self.protein_library[uid] = metadata

            print(f"Library built with {len(self.protein_library)} proteins.")

        except Exception as e:
            print(f"Error fetching data: {e}")

    def export_library_to_excel(engine, filename="SPARTA_Metadata_Sanity_Check.xlsx"):
        """
        Converts the Protein Library into a DataFrame and exports to Excel.
        
        Processing:
        - Flattens PDB ID lists into comma-separated strings.
        - Simplifies MS Comments for spreadsheet readability.
        """
        if not engine.protein_library:
            print("Error: Library is empty. Fetch data before exporting.")
            return

        # 1. Convert the dictionary of dataclasses into a list of dictionaries
        raw_data = []
        for uid, metadata in engine.protein_library.items():
            # Convert dataclass to dict
            entry_dict = asdict(metadata)
            
            # 2. Format lists for Excel (Excel cells don't handle Python lists well)
            entry_dict['pdb_ids'] = ", ".join(entry_dict['pdb_ids'])
                
            raw_data.append(entry_dict)

        # 4. Create DataFrame
        df = pd.DataFrame(raw_data)

        # 5. Reorder columns for easier "Sanity Checking"
        cols = [
            'uniprot_id', 'entry_name', 'protein_name', 
            'mono_mass', 'mono_mass_no_met', 
            'avg_mass', 'avg_mass_no_met',
            'hierarchy_priority', 'pdb_ids', 'alphafold_id'
        ]
        df = df[cols]

        # 6. Export
        df.to_excel(filename, index=False)
        print(f"✅ Sanity check file generated: {filename}")


# # Example Usage:
# engine = MetadataEngine()
# engine.fetch_proteome(10090) # Mouse
# engine.export_library_to_excel()
# hits = engine.search_by_mass(257620.0, ppm_tolerance=20.0) # Search for Ubiquitin
# print(hits)
