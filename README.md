# GeneWeaver – DNA Sequence Analysis System

GeneWeaver is an educational Python web application that validates and analyzes DNA sequences.

## Features
- DNA sequence validation (A, T, G, C)
- Nucleotide counts and sequence length
- GC and AT content percentages
- Reverse complement
- DNA-to-RNA transcription
- Basic codon translation to one-letter amino-acid codes
- FASTA and plain-text file upload

## Requirements
- Python 3.10 or newer
- pip

## Run on Windows (VS Code PowerShell)

Open the `GeneWeaver` folder in VS Code. In Terminal → New Terminal, run:

```powershell
py --version
py -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
pip install -r requirements.txt
python app.py
```

If PowerShell blocks virtual-environment activation, run this once in that terminal:

```powershell
Set-ExecutionPolicy -Scope Process -ExecutionPolicy Bypass
.\.venv\Scripts\Activate.ps1
```

Then open this address in Chrome:

http://127.0.0.1:5000

Stop the server with `Ctrl + C` in the terminal.

## Test example
Paste this sequence into the form:

`ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG`

Click **Analyze DNA**.

## Project structure

```text
GeneWeaver/
├── app.py
├── analyzer.py
├── requirements.txt
├── README.md
├── templates/
│   └── index.html
└── static/
    └── style.css
```

## Notes
This is a learning/demo tool. Translation starts at the first base, reads complete codons, and stops at the first stop codon. It does not perform clinical or medical interpretation.
