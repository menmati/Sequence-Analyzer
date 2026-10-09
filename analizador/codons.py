"""
Módulo para el análisis del contenido de GC, codones de inicio/parada y marcos abiertos de lectura (ORFs).
"""

from dataclasses import dataclass
from typing import Dict, List, Tuple


START_CODONS = {"ATG", "AUG"}
STOP_CODONS = {"TAA", "TAG", "TGA", "UAA", "UAG", "UGA"}


@dataclass
class CodonMatch:
    """Información de un codón localizado en la secuencia."""
    codon: str
    kind: str           # 'inicio' o 'parada'
    position: int       # Posición base 1 en la secuencia
    frame: int          # Marco de lectura (1, 2 o 3)


@dataclass
class ORF:
    """Marco abierto de lectura (Open Reading Frame)."""
    frame: int
    start_pos: int      # Posición base 1 del codón de inicio
    stop_pos: int       # Posición base 1 del final del codón de parada
    length_nt: int      # Longitud en nucleótidos
    length_aa: int      # Longitud en aminoácidos (excluyendo o incluyendo stop)
    sequence_nt: str    # Secuencia de nucleótidos del ORF


def calculate_gc_content(sequence: str) -> Dict[str, float]:
    """
    Calcula el porcentaje de contenido GC y frecuencias de nucleótidos.

    Returns:
        Diccionario con porcentaje de GC, y conteo/porcentajes de A, C, G, T.
    """
    seq = sequence.upper().replace("U", "T")
    total_len = len(seq)
    if total_len == 0:
        return {
            "gc_percent": 0.0,
            "total_bases": 0,
            "counts": {"A": 0, "C": 0, "G": 0, "T": 0, "otros": 0},
            "frequencies": {"A": 0.0, "C": 0.0, "G": 0.0, "T": 0.0}
        }

    counts = {
        "A": seq.count("A"),
        "C": seq.count("C"),
        "G": seq.count("G"),
        "T": seq.count("T"),
    }
    standard_bases = sum(counts.values())
    counts["otros"] = total_len - standard_bases

    gc_count = counts["G"] + counts["C"]
    # Si hay bases ambiguas (N, etc.), calculamos sobre bases conocidas o total
    gc_percent = (gc_count / standard_bases * 100.0) if standard_bases > 0 else 0.0

    frequencies = {
        base: (count / standard_bases * 100.0) if standard_bases > 0 else 0.0
        for base, count in counts.items() if base != "otros"
    }

    return {
        "gc_percent": round(gc_percent, 2),
        "total_bases": total_len,
        "standard_bases": standard_bases,
        "counts": counts,
        "frequencies": {k: round(v, 2) for k, v in frequencies.items()}
    }


def find_codons(sequence: str, frame: int = 0) -> Tuple[List[CodonMatch], List[CodonMatch]]:
    """
    Encuentra codones de inicio y parada en la secuencia.

    Args:
        sequence: Cadena de nucleótidos.
        frame: 0 para buscar en los 3 marcos (1, 2, 3), o especificar 1, 2 o 3.

    Returns:
        Tupla (start_codons, stop_codons) con listas de CodonMatch.
    """
    seq = sequence.upper().replace("U", "T")
    start_matches: List[CodonMatch] = []
    stop_matches: List[CodonMatch] = []

    frames_to_check = [1, 2, 3] if frame == 0 else [frame]

    for f in frames_to_check:
        start_idx = f - 1
        for i in range(start_idx, len(seq) - 2, 3):
            codon = seq[i:i+3]
            pos_1based = i + 1
            if codon in {"ATG"}:
                start_matches.append(CodonMatch(codon=codon, kind="inicio", position=pos_1based, frame=f))
            elif codon in {"TAA", "TAG", "TGA"}:
                stop_matches.append(CodonMatch(codon=codon, kind="parada", position=pos_1based, frame=f))

    # Ordenar por posición
    start_matches.sort(key=lambda x: x.position)
    stop_matches.sort(key=lambda x: x.position)

    return start_matches, stop_matches


def find_orfs(sequence: str, min_length_nt: int = 30, frame: int = 0) -> List[ORF]:
    """
    Identifica Open Reading Frames (ORFs) entre un codón de inicio (ATG)
    y el siguiente codón de parada en el mismo marco de lectura.

    Args:
        sequence: Secuencia de nucleótidos.
        min_length_nt: Longitud mínima del ORF en nucleótidos (por defecto 30 nt = 10 aa).
        frame: 0 para los 3 marcos, o 1, 2, 3.

    Returns:
        Lista de instancias de ORF ordenadas por posición de inicio.
    """
    seq = sequence.upper().replace("U", "T")
    orfs: List[ORF] = []
    frames_to_check = [1, 2, 3] if frame == 0 else [frame]

    for f in frames_to_check:
        start_idx = f - 1
        i = start_idx
        while i <= len(seq) - 3:
            codon = seq[i:i+3]
            if codon == "ATG":
                start_pos = i
                # Buscar el siguiente codón de parada en este marco
                stop_found = False
                for j in range(i + 3, len(seq) - 2, 3):
                    stop_codon = seq[j:j+3]
                    if stop_codon in {"TAA", "TAG", "TGA"}:
                        orf_seq = seq[start_pos:j+3]
                        if len(orf_seq) >= min_length_nt:
                            orfs.append(ORF(
                                frame=f,
                                start_pos=start_pos + 1,
                                stop_pos=j + 3,
                                length_nt=len(orf_seq),
                                length_aa=len(orf_seq) // 3,
                                sequence_nt=orf_seq
                            ))
                        # Continuar la búsqueda desde después de este inicio
                        stop_found = True
                        break
            i += 3

    orfs.sort(key=lambda x: (x.start_pos, -x.length_nt))
    return orfs
