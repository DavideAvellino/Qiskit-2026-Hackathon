# Grover × DNA: screening di mutazioni con Qiskit

Tool che cerca una mutazione genetica nel DNA di un paziente con l'algoritmo di Grover.
Esempio reale: anemia falciforme, gene HBB (NCBI NM_000518.5).

## Contenuto

```
codice/
  screening_grover.py         il tool, tutto in un file
  grover_dna_spiegato.ipynb   notebook passo passo (oracolo, diffusore, versione con le lettere)
  requirements.txt            librerie necessarie
presentazione/
  presentazione.pdf           slide del pitch (9 slide, 5 minuti)
  presentazione.html          le stesse slide, apribili nel browser
  note_pitch.md               cosa dire per ogni slide + domande probabili della giuria
spiegazioni.md                riassunto dei concetti da studiare
```

## Come far girare il tool

```bash
cd codice
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
python screening_grover.py
```
