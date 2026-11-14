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
- **Sequence Transformations**:
  - 5' → 3' Reverse complement.
  - Transcription from DNA to mRNA ($T \rightarrow U$).
  - Full codon translation into IUPAC amino acid single-letter codes with premature stop codon detection.
- **Open Reading Frame (ORF) Detection**: Automated discovery of in-frame protein coding segments starting with `ATG` across forward reading frames.
- **Multi-FASTA Support**: Robust parsing of single or multiple FASTA header records.
- **Export Options**: Download comprehensive reports as structured **JSON** or formatted **TXT**.
- **Dual Interfaces**: Full-featured **Flask Web Application** and a scriptable **CLI Tool** (`cli.py`).
- **Developer REST API**: Programmatic endpoints (`/api/analyze`, `/api/health`) for pipeline automation.

---

## 📁 Project Structure

```text
GeneWeaver/
├── app.py                  # Flask web server, REST API & download routes
├── analyzer.py             # Core bioinformatics analysis algorithms & calculations
├── cli.py                  # Command-line interface with argparse
├── requirements.txt        # Python dependency manifest
├── README.md               # Project documentation
├── sample.fasta            # Sample FASTA test record
├── tests/
│   └── test_analyzer.py    # Unit tests suite (17 test cases)
├── templates/
│   └── index.html          # Interactive Jinja2 template with gel simulator & charts
└── static/
    └── style.css           # Modern responsive stylesheet with agarose gel theme
```

---

## 🛠️ Setup & Execution

### 1. Environment Setup (Windows PowerShell)

```powershell
# Open the project folder in terminal
cd d:\GeneWeaver

# Create and activate virtual environment
python -m venv .venv
.\.venv\Scripts\Activate.ps1

# Install requirements
pip install -r requirements.txt
```

### 2. Run the Web Application

```powershell
python app.py
```

Open your browser at: **`http://127.0.0.1:5000`**

---

## 💻 Command Line Interface (CLI)

GeneWeaver includes a command-line tool (`cli.py`):

```powershell
# 1. Analyze an inline sequence
python cli.py -s ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG

# 2. Analyze a FASTA file and output as JSON
python cli.py -f sample.fasta --json

# 3. Detect Open Reading Frames exclusively
python cli.py -f sample.fasta --orfs

# 4. Show Restriction Enzyme Digestion Profile
python cli.py -f sample.fasta --enzymes

# 5. Save analysis report to a file
python cli.py -f sample.fasta -o report.txt
```

---

## 🔌 REST API Documentation

### Health Check
`GET /api/health`
```json
{
  "status": "healthy",
  "service": "GeneWeaver DNA Analysis API",
  "version": "2.1.0"
}
```

### Sequence Analysis
`POST /api/analyze`
```bash
curl -X POST http://127.0.0.1:5000/api/analyze \
  -H "Content-Type: application/json" \
  -d '{"sequence": "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"}'
```

---

## 🧪 Running Unit Tests

Run the complete test suite (17/17 passing):

```powershell
python -m unittest discover tests
```

---

## 📜 Disclaimer
This software is developed as an educational demonstration project. It is intended for computational learning, academic research simulations, and bioinformatics training.
