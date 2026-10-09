"""
Módulo para la traducción de secuencias de ADN/ARN a secuencias proteicas.
"""

from typing import Dict, List, Optional
from .fasta import reverse_complement


# Código genético estándar (Tabla NCBI 1)
GENETIC_CODE: Dict[str, str] = {
    # Fenilalanina
    "TTT": "F", "TTC": "F",
    # Leucina
    "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    # Isoleucina
    "ATT": "I", "ATC": "I", "ATA": "I",
    # Metionina (Inicio)
    "ATG": "M",
    # Valina
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    # Serina
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "AGT": "S", "AGC": "S",
    # Prolina
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    # Treonina
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    # Alanina
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    # Tirosina
    "TAT": "Y", "TAC": "Y",
    # Codones de parada
    "TAA": "*", "TAG": "*", "TGA": "*",
    # Histidina
    "CAT": "H", "CAC": "H",
    # Glutamina
    "CAA": "Q", "CAG": "Q",
    # Asparagina
    "AAT": "N", "AAC": "N",
    # Lisina
    "AAA": "K", "AAG": "K",
    # Ácido Aspártico
    "GAT": "D", "GAC": "D",
    # Ácido Glutámico
    "GAA": "E", "GAG": "E",
    # Cisteína
    "TGT": "C", "TGC": "C",
    # Triptófano
    "TGG": "W",
    # Arginina
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
    # Glicina
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

# Nombres completos de aminoácidos
AMINO_ACID_NAMES: Dict[str, str] = {
    "A": "Alanina (Ala)",
    "R": "Arginina (Arg)",
    "N": "Asparagina (Asn)",
    "D": "Ácido Aspártico (Asp)",
    "C": "Cisteína (Cys)",
    "E": "Ácido Glutámico (Glu)",
    "Q": "Glutamina (Gln)",
    "G": "Glicina (Gly)",
    "H": "Histidina (His)",
    "I": "Isoleucina (Ile)",
    "L": "Leucina (Leu)",
    "K": "Lisina (Lys)",
    "M": "Metionina (Met / Inicio)",
    "F": "Fenilalanina (Phe)",
    "P": "Prolina (Pro)",
    "S": "Serina (Ser)",
    "T": "Treonina (Thr)",
    "W": "Triptófano (Trp)",
    "Y": "Tirosina (Tyr)",
    "V": "Valina (Val)",
    "*": "Codón de Parada (Stop)",
    "X": "Desconocido (Xaa)",
}


def translate_codon(codon: str) -> str:
    """Traduce un codón de 3 nucleótidos a su aminoácido correspondiente."""
    codon_clean = codon.upper().replace("U", "T")
    return GENETIC_CODE.get(codon_clean, "X")


def translate_dna(
    sequence: str,
    frame: int = 1,
    stop_symbol: str = "*",
    to_first_stop: bool = False
) -> str:
    """
    Traduce una secuencia nucleotídica a proteína en un marco de lectura dado.

    Args:
        sequence: Secuencia de ADN/ARN.
        frame: Marco de lectura (1, 2, o 3).
        stop_symbol: Símbolo para codones de parada (por defecto '*').
        to_first_stop: Si es True, detiene la traducción al encontrar el primer codón de parada.

    Returns:
        Cadena de aminoácidos resultante.
    """
    if frame not in {1, 2, 3}:
        raise ValueError(f"El marco de lectura debe ser 1, 2 o 3. Se recibió: {frame}")

    seq = sequence.upper().replace("U", "T")
    start_idx = frame - 1
    protein_chars: List[str] = []

    for i in range(start_idx, len(seq) - 2, 3):
        codon = seq[i:i+3]
        aa = translate_codon(codon)
        if aa == "*":
            if to_first_stop:
                protein_chars.append(stop_symbol)
                break
            protein_chars.append(stop_symbol)
        else:
            protein_chars.append(aa)

    return "".join(protein_chars)


def translate_all_frames(sequence: str, include_reverse: bool = False) -> Dict[str, str]:
    """
    Traduce la secuencia en múltiples marcos de lectura.

    Args:
        sequence: Secuencia nucleotídica.
        include_reverse: Si es True, también traduce marcos -1, -2, -3 de la hebra complementaria reversa.

    Returns:
        Diccionario con la clave del marco (+1, +2, +3, -1, -2, -3) y la proteína correspondiente.
    """
    results: Dict[str, str] = {}
    for f in [1, 2, 3]:
        results[f"+{f}"] = translate_dna(sequence, frame=f)

    if include_reverse:
        rev_seq = reverse_complement(sequence)
        for f in [1, 2, 3]:
            results[f"-{f}"] = translate_dna(rev_seq, frame=f)

    return results
