# GeneWeaver – DNA Sequence Analysis & Bioinformatics Workbench

GeneWeaver is a modular Python web application and command-line bioinformatics suite designed to analyze, validate, and explore DNA sequences.

---

## 🚀 Key Features

- **Sequence Validation**: Strict IUPAC base verification (`A`, `T`, `G`, `C`) with whitespace normalization.
- **Biochemical Properties**:
  - Nucleotide composition counts and percentage ratios.
  - GC / AT content percentages.
  - **GC Skew** calculation: $(G - C) / (G + C)$ across the sequence.
  - **Melting Temperature ($T_m$)**: Wallace empirical rule for short oligos and salt-adjusted formula for extended sequences.
  - **Molecular Weight ($M_w$)**: ssDNA and dsDNA mass calculation in g/mol.
  - **Protein Isoelectric Point ($pI$) & Net Charge**: Bisection search for zero-net-charge pH and charge estimation at physiological pH 7.4.
- **Restriction Endonuclease Digestion & Gel Simulation**:
  - Recognition site scanning for standard cloning enzymes (`EcoRI`, `BamHI`, `HindIII`, `NotI`, `PstI`, `SmaI`, `XhoI`, `TaqI`, `HaeIII`, `AluI`).
  - Cut coordinates and digestion fragment length calculation.
  - **In Silico Agarose Gel Electrophoresis**: Visual lane simulation comparing sample digestion bands against a 1 kb DNA ladder.
- **Pairwise Alignment & SNP Variant Detection**:
  - Detects point mutations between query/mutant and reference sequences.
  - Transition vs. Transversion classification with genomic $T_s / T_v$ quality ratio.
  - Functional impact prediction: Synonymous (Silent), Missense, and Nonsense (Stop gained).
- **PCR Primer Design & Reverse-Translation**:
  - Automated flanking forward and reverse PCR primer design with $T_m$ compatibility checking.
  - Reverse-translation of peptide sequences into consensus degenerate DNA using IUPAC ambiguity codes.
- **Multi-FASTA Batch Analysis & Genetic Distance Matrix**:
  - High-throughput comparative analysis across multi-sequence FASTA files.
  - Pairwise $N \times N$ genetic distance matrix (Hamming distance and percent identity).
  - Automated detection of closest homologous pairs and most divergent outliers.
- **Dual Interfaces & REST API**:
  - Interactive **Flask Web Application** with tabbed navigation.
  - Scriptable **CLI Tool** (`cli.py`).
  - Programmatic **REST API** (`/api/analyze`, `/api/compare`, `/api/primers`, `/api/reverse-translate`, `/api/batch-analyze`, `/api/health`).

---

## 📁 Project Structure

```text
GeneWeaver/
├── app.py                  # Flask web server, REST API & download routes
├── analyzer.py             # Core bioinformatics analysis algorithms & calculations
├── batch_analyzer.py       # Multi-FASTA comparative metrics & distance matrix
├── primers.py              # PCR primer design & reverse-translation engine
├── cli.py                  # Command-line interface with argparse
├── requirements.txt        # Python dependency manifest
├── README.md               # Project documentation
├── sample.fasta            # Sample single-record FASTA
├── sample_multi.fasta      # Sample multi-record comparative FASTA
├── tests/
│   ├── test_analyzer.py        # Core algorithm unit tests (20 tests)
│   ├── test_primers.py         # Primer design & reverse-translation unit tests (5 tests)
│   └── test_batch_analyzer.py  # Batch analysis & distance matrix unit tests (5 tests)
├── templates/
│   └── index.html          # Interactive Jinja2 template with gel simulator & tabs
└── static/
    └── style.css           # Modern responsive stylesheet
```

---

## 🛠️ Setup & Execution

### 1. Environment Setup (Windows PowerShell)

```powershell
cd d:\GeneWeaver
.\.venv\Scripts\Activate.ps1
pip install -r requirements.txt
```

### 2. Run the Web Application

```powershell
python app.py
```
Open your browser at: **`http://127.0.0.1:5000`**

---

## 💻 Command Line Interface (CLI)

```powershell
# 1. Analyze single sequence
python cli.py -s ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG

# 2. Compare sequence against reference (SNP variant detection)
python cli.py -s ATGGCCATTGTA --compare ATGGCCATGGTA

# 3. Flanking PCR primer design
python cli.py -f sample.fasta --primers

# 4. Restriction enzyme digestion profile
python cli.py -f sample.fasta --enzymes

# 5. Multi-FASTA batch analysis & genetic distance matrix
python cli.py -f sample_multi.fasta --batch
```

---

## 🧪 Running Unit Tests

Run the complete test suite (30/30 tests passing):

```powershell
python -m unittest discover tests
```

---

## 📜 Disclaimer
This software is developed as an educational demonstration project for computational learning, academic research simulations, and bioinformatics training.
