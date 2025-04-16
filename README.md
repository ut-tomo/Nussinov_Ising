# RNA Secondary Structure Prediction

This project demonstrates two methods for predicting RNA secondary structures:

1. **Quantum-inspired optimization** using the [Fixstars Amplify](https://amplify.fixstars.com/) SDK
2. **Classical dynamic programming** using the Nussinov algorithm

---

## 🔬 Overview

- Randomly generates an RNA sequence of length `n`
- Computes valid base pairings based on biological pairing rules (including wobble pairs)
- Constructs an optimization problem with constraints:
  - Each base pairs at most once
  - No crossing pairs
  - Minimum loop length of 3 bases
- Solves the problem using:
  - Amplify (QUBO formulation)
  - Nussinov (dynamic programming)
- Compares predicted secondary structures and scores

---

## 🧬 RNA Pairing Rules

The valid base pairs include:

- **A–U** and **U–A**
- **G–C** and **C–G**
- **G–U** and **U–G** (wobble pairs)

---

## ⚙️ Requirements

- Python 3.9~
- `numpy`
- `amplify` (Fixstars Amplify SDK)

You can install dependencies via:

```bash
pip install numpy amplify
