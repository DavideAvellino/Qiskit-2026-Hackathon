"""
Screening di una mutazione genetica con l'algoritmo di Grover.

Esempio: anemia falciforme, gene HBB (fonte NCBI NM_000518.5).
Uso: python screening_grover.py

"""

import math

from qiskit import QuantumCircuit
from qiskit.primitives import StatevectorSampler


# ---------------------------------------------------------------------------
# 1. I DATI
# ---------------------------------------------------------------------------

# Gene HBB di una persona sana (444 lettere)
DNA_SANO = (
    "ATGGTGCATCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAG"
    "TTGGTGGTGAGGCCCTGGGCAGGCTGCTGGTGGTCTACCCTTGGACCCAGAGGTTCTTTGAGTCCTTTGG"
    "GGATCTGTCCACTCCTGATGCTGTTATGGGCAACCCTAAGGTGAAGGCTCATGGCAAGAAAGTGCTCGGT"
    "GCCTTTAGTGATGGCCTGGCTCACCTGGACAACCTCAAGGGCACCTTTGCCACACTGAGTGAGCTGCACT"
    "GTGACAAGCTGCACGTGGATCCTGAGAACTTCAGGCTCCTGGGCAACGTGCTGGTCTGTGTGCTGGCCCA"
    "TCACTTTGGCAAAGAATTCACCCCACCAGTGCAGGCTGCCTATCAGAAAGTGGTGGCTGGTGTGGCTAAT"
    "GCCCTGGCCCACAAGTATCACTAA"
)

# Anemia falciforme: una sola lettera cambia (indice 19, da A a T), quindi GAG diventa GTG
DNA_MALATO = DNA_SANO[:19] + "T" + DNA_SANO[20:]

# La mutazione da cercare: la lettera cambiata con un po' di contesto attorno
MUTAZIONE = "CCTGTGGAG"


# ---------------------------------------------------------------------------
# 2. I PEZZI DEL CIRCUITO
# ---------------------------------------------------------------------------

def z_multicontrollata(circuito):
    # Cambia il segno solo dello stato |11...1> (tutti i qubit a 1)
    ultimo = circuito.num_qubits - 1
    circuito.h(ultimo)
    circuito.mcx(list(range(ultimo)), ultimo)
    circuito.h(ultimo)


def crea_oracolo(numero_qubit, posizioni_giuste):
    # Mette il segno meno sugli indici dove comincia la mutazione
    oracolo = QuantumCircuit(numero_qubit)
    for posizione in posizioni_giuste:
        # Qubit che valgono 0 in questa posizione: con X li portiamo a 1
        qubit_a_zero = [q for q in range(numero_qubit) if (posizione >> q) & 1 == 0]
        if qubit_a_zero:
            oracolo.x(qubit_a_zero)
        z_multicontrollata(oracolo)  # segno meno su |11...1>
        if qubit_a_zero:
            oracolo.x(qubit_a_zero)  # rimettiamo tutto com'era
    return oracolo


def crea_diffusore(numero_qubit):
    # "Inversione attorno alla media": fa crescere gli stati con il segno meno
    diffusore = QuantumCircuit(numero_qubit)
    diffusore.h(range(numero_qubit))
    diffusore.x(range(numero_qubit))
    z_multicontrollata(diffusore)
    diffusore.x(range(numero_qubit))
    diffusore.h(range(numero_qubit))
    return diffusore


# ---------------------------------------------------------------------------
# 3. LA RICERCA CON GROVER
# ---------------------------------------------------------------------------

def cerca_mutazione(dna, mutazione):
    n = len(dna)
    m = len(mutazione)
    numero_posizioni = n - m + 1  # indici possibili: da 0 a n - m

    # Posizioni che l'oracolo deve segnare (le calcola Python)
    posizioni_giuste = [i for i in range(numero_posizioni) if dna[i : i + m] == mutazione]

    # Qubit: quelli per scrivere l'indice + 2 in più (Grover funziona meglio)
    numero_qubit = math.ceil(math.log2(numero_posizioni)) + 2
    numero_stati = 2 ** numero_qubit

    # Ripetizioni ottimali: pi/4 * radice(N / k). Se la mutazione non c'è usiamo k = 1 (evita la divisione per zero)
    k = max(1, len(posizioni_giuste))
    ripetizioni = math.floor(math.pi / 4 * math.sqrt(numero_stati / k))

    # Circuito: sovrapposizione, poi (oracolo + diffusore) ripetuti, poi misura
    circuito = QuantumCircuit(numero_qubit)
    circuito.h(range(numero_qubit))
    oracolo = crea_oracolo(numero_qubit, posizioni_giuste)
    diffusore = crea_diffusore(numero_qubit)
    for _ in range(ripetizioni):
        circuito.compose(oracolo, inplace=True)
        circuito.compose(diffusore, inplace=True)
    circuito.measure_all()

    # Eseguiamo 1024 volte sul simulatore e prendiamo l'indice uscito più spesso
    conteggi = StatevectorSampler().run([circuito], shots=1024).result()[0].data.meas.get_counts()
    bit_piu_frequenti = max(conteggi, key=conteggi.get)
    indice = int(bit_piu_frequenti, 2)
    frequenza = conteggi[bit_piu_frequenti] / 1024

    # Verifica classica: in quell'indice c'è davvero la mutazione?
    trovata = dna[indice : indice + m] == mutazione

    print(f"DNA analizzato:     {n} lettere")
    print(f"Mutazione cercata:  {mutazione}")
    print(f"Qubit usati:        {numero_qubit}")
    print(f"Ripetizioni Grover: {ripetizioni}")
    if trovata:
        print(f"Indice proposto:    {indice} (uscito nel {frequenza:.1%} delle misure, verifica OK)")
        print(f">>> MUTAZIONE TROVATA all'indice {indice}: consigliati ulteriori controlli")
    else:
        print(f"Indice proposto:    nessuno si distingue (il più frequente solo nel {frequenza:.1%})")
        print(">>> Mutazione NON presente")


# ---------------------------------------------------------------------------
# 4. LA DEMO
# ---------------------------------------------------------------------------

print("Screening anemia falciforme - gene HBB\n")
print("=== Paziente 1 ===")
cerca_mutazione(DNA_SANO, MUTAZIONE)
print("\n=== Paziente 2 ===")
cerca_mutazione(DNA_MALATO, MUTAZIONE)
