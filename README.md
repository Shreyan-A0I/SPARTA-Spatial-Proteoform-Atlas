# SPARTA: Spatial Proteoform Atlas

**Multi-modal discovery platform integrating MALDI Mass Spectrometry Imaging (MSI) with structural bioinformatics to map Activity Zones within tissues.**

## 🎯 Project Overview

The Spatial Proteoform Atlas moves beyond simple protein identification to reveal the functional landscape of the proteome in its native spatial context. By correlating spatially resolved mass shifts (indicating Post-Translational Modifications) with experimental 3D structures from the PDB, the platform maps "Activity Zones" within tissues.

## 📋 Technical Workflow

### Phase A: The Metadata Engine (The Blueprint)
Constructs a theoretical mass library with structural cross-references.

**Mathematical Foundations:**
- **Mass Calculation**: $m_{th} = \sum (\text{Residue Masses}) + 18.01056 \text{ Da}$
- **Methionine Logic**: Stores two potential masses: $m_{th}$ and $(m_{th} - 131.0405 \text{ Da})$
- **Hierarchy Logic**: PDB (X-ray/Cryo-EM) > AlphaFold (Predicted)

### Phase B: The Spectral Pipeline (The Evidence)
Denoises raw MSI data and normalizes across the tissue slice.

**Mathematical Foundations:**
- **TIC Normalization**: $I_{norm} = \frac{I_i}{\sum_{j=1}^{n} I_j}$
- **Peak Picking**: $S/N = \frac{I_{peak} - \mu_{noise}}{\sigma_{noise}} \geq 3$
- **Cosine Similarity**: $\text{similarity} = \cos(\theta) = \frac{\mathbf{A} \cdot \mathbf{B}}{\|\mathbf{A}\| \|\mathbf{B}\|}$

### Phase C: The Discovery Matcher (The Bridge)
Identifies functional mass shifts (PTMs) with statistical confidence.

**Mathematical Foundations:**
- **Error Measurement**: $E_{ppm} = \frac{|m_{obs} - m_{th}|}{m_{th}} \times 10^6 \leq 50 \text{ ppm}$
- **Delta Matcher**: $m_{obs} \approx (m_{th} + \Delta_{PTM})$
- **Pearson Correlation**: $r = \frac{\sum (X_i - \bar{X})(Y_i - \bar{Y})}{\sqrt{\sum (X_i - \bar{X})^2 \sum (Y_i - \bar{Y})^2}}$
- **Confidence Threshold**: $r > 0.6$ suggests high-probability proteoform

### Phase D: Structural Projection (The Visualization)
Maps discovered modifications to 3D physical models.

**Mathematical Foundations:**
- **SASA Calculation**: Determines if a residue is physically reachable
- **Accessibility Filter**: Residue is "Likely Target" if $\text{SASA} > 20 \text{ \AA}^2$
- **Confidence Scoring**: $C = w_1(E_{ppm}) + w_2(r) + w_3(\text{SASA})$

## 📁 Repository Structure

```
SPARTA/
├── data/                   # Raw MSI data (Keep out of Git)
│   ├── mouse.ibd
│   └── mouse.imzML
├── src/                    # All source code
│   ├── backend/
│   │   ├── metadata_engine.py  # Task T1: UniProt/PDB scraper
│   │   ├── spectral_proc.py    # Task T2: MSI processing logic
│   │   └── discovery_logic.py  # Task T3/T4: Delta-matching & Stats
│   ├── frontend/
│   │   └── dashboard.py        # Task T6: Streamlit app
│   └── utils/
│       └── sasa_engine.py      # Task T5: Structural calculations
├── notebooks/              # For R&D and testing algorithms
├── requirements.txt        # Tech stack dependencies
└── README.md               # Project documentation
```

## 🚀 Installation

1. **Clone the repository:**
   ```bash
   git clone <repository-url>
   cd SPARTA
   ```

2. **Create a virtual environment:**
   ```bash
   python -m venv venv
   source venv/bin/activate  # On Windows: venv\Scripts\activate
   ```

3. **Install dependencies:**
   ```bash
   pip install -r requirements.txt
   ```

## 💻 Usage

### Running the Dashboard

```bash
streamlit run src/frontend/dashboard.py
```

### Using Individual Modules

```python
from src.backend.metadata_engine import MetadataEngine
from src.backend.spectral_proc import SpectralProcessor
from src.backend.discovery_logic import DiscoveryMatcher
from src.utils.sasa_engine import SASAEngine

# Initialize engines
metadata = MetadataEngine()
spectral = SpectralProcessor()
discovery = DiscoveryMatcher()
sasa = SASAEngine()

# Your workflow here...
```

## 📊 Task Roadmap

| Task ID | Task Description | Math/Logic Focus | Hours |
|---------|------------------|------------------|-------|
| T1 | UniProt Scraper | Hierarchy & Met-loss logic | 15h |
| T2 | imzML Loader | TIC Normalization & Peak Picking | 25h |
| T3 | PTM Logic | ppm Error & Delta-Mass matching | 20h |
| T4 | Spatial Stats | Pearson $r$ & Co-localization | 15h |
| T5 | SASA Engine | 3D coordinate mapping & Surface Area | 20h |
| T6 | UI/UX Dashboard | Interactive Heatmaps & 3D rendering | 40h |

**Total Estimated Time: ~135 Hours**

## 🛠️ Tech Stack

- **Data Handling**: Pandas, NumPy, H5Py
- **Signal Processing**: pyOpenMS, pyimzML
- **Bioinformatics**: Biopython (PDB parsing), FreeSASA
- **Frontend**: Streamlit, Plotly (spatial heatmaps), Stmol (NGLview integration)
- **Statistics**: SciPy, scikit-learn

## ⚠️ Risks & Mitigations

### Isobaric Overlap
**Risk**: 80 Da can be Phosphorylation or Sulfation.

**Mitigation**: The UI flags "Isobaric Ambiguity" when mass differences are $< 0.01 \text{ Da}$.

### Fragmentary PDBs
**Risk**: PDB files often represent only a domain of a protein.

**Mitigation**: Implement Sequence Alignment (using Biopython's Bio.Align). If PDB coverage $< 70\%$, the app automatically triggers an AlphaFold DB fallback for the full-length structure.

## 📝 Development Status

This project is currently in **Phase 1: Foundation Setup**. All core modules have been scaffolded with placeholder implementations. Next steps:

1. ✅ Project structure created
2. ✅ Core modules scaffolded
3. ⏳ Task T1: Implement UniProt API integration
4. ⏳ Task T2: Implement imzML loading and processing
5. ⏳ Task T3: Implement PTM matching logic
6. ⏳ Task T4: Implement spatial statistics
7. ⏳ Task T5: Implement SASA calculations
8. ⏳ Task T6: Complete dashboard implementation

## 🤝 Contributing

Contributions are welcome! Please follow these guidelines:

1. Create a feature branch
2. Implement your changes with appropriate tests
3. Update documentation as needed
4. Submit a pull request

## 📄 License

See LICENSE file for details.

## 🙏 Acknowledgments

- UniProt Consortium for protein data
- Protein Data Bank (PDB) for structural data
- AlphaFold Database for predicted structures
- OpenMS community for mass spectrometry tools

