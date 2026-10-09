"""
Módulo para la detección y clasificación de mutaciones puntuales (SNPs y sustituciones).
"""

from dataclasses import dataclass
from typing import List, Optional, Tuple
from .translation import translate_codon, AMINO_ACID_NAMES


PURINES = {"A", "G"}
PYRIMIDINES = {"C", "T"}


@dataclass
class PointMutation:
    """Representa una mutación puntual detectada entre dos secuencias."""
    position: int                  # Posición en base 1
    ref_base: str                  # Base en la secuencia de referencia
    alt_base: str                  # Base en la secuencia mutada
    type: str                      # 'Transición', 'Transversión', o 'Indel/Gaps'
    codon_pos: Optional[int]       # Posición del codón (1-based en la proteína)
    ref_codon: Optional[str]       # Codón de referencia
    alt_codon: Optional[str]       # Codón mutado
    ref_aa: Optional[str]          # Aminoácido referencia
    alt_aa: Optional[str]          # Aminoácido mutado
    effect: Optional[str]          # 'Sinónima (Silenciosa)', 'No sinónima (Missense)', 'Sin sentido (Nonsense)', etc.
    description: str               # Descripción textual explicativa


def classify_nucleotide_substitution(ref: str, alt: str) -> str:
    """Clasifica una sustitución nucleotídica en Transición o Transversión."""
    ref = ref.upper().replace("U", "T")
    alt = alt.upper().replace("U", "T")

    if (ref in PURINES and alt in PURINES) or (ref in PYRIMIDINES and alt in PYRIMIDINES):
        return "Transición"
    elif (ref in PURINES and alt in PYRIMIDINES) or (ref in PYRIMIDINES and alt in PURINES):
        return "Transversión"
    return "Indel / Inespecífica"


def classify_protein_effect(ref_aa: str, alt_aa: str) -> str:
    """Determina el efecto funcional de la mutación a nivel proteico."""
    if ref_aa == alt_aa:
        return "Sinónima (Silenciosa)"
    elif ref_aa != "*" and alt_aa == "*":
        return "Sin sentido (Nonsense / Parada prematura)"
    elif ref_aa == "*" and alt_aa != "*":
        return "Pérdida de parada (Readthrough)"
    elif ref_aa == "M" and alt_aa != "M":
        return "Pérdida de codón de inicio"
    else:
        return "No sinónima (Missense / Cambio de sentido)"


def get_codon_at_position(seq: str, pos_0based: int, frame: int = 1) -> Tuple[Optional[str], Optional[int], int]:
    """
    Obtiene el triplete (codón) que contiene una posición dada según el marco de lectura.

    Returns:
        Tupla (triplete, posicion_del_codon_en_proteina, offset_en_el_codon 0..2)
    """
    start_frame = frame - 1
    if pos_0based < start_frame:
        return None, None, -1

    rel_pos = pos_0based - start_frame
    codon_idx = rel_pos // 3
    codon_offset = rel_pos % 3
    codon_start = start_frame + (codon_idx * 3)

    if codon_start + 3 <= len(seq):
        return seq[codon_start:codon_start + 3], codon_idx + 1, codon_offset
    return None, None, -1


def identify_point_mutations(
    ref_sequence: str,
    query_sequence: str,
    frame: int = 1
) -> List[PointMutation]:
    """
    Compara la secuencia problema (query) contra la secuencia de referencia
    e identifica todas las mutaciones puntuales.

    Args:
        ref_sequence: Secuencia nucleotídica de referencia.
        query_sequence: Secuencia nucleotídica comparada.
        frame: Marco de lectura para determinar el impacto en la traducción (1, 2 o 3).

    Returns:
        Lista de mutaciones puntuales detectadas.
    """
    ref = ref_sequence.upper().replace("U", "T")
    query = query_sequence.upper().replace("U", "T")

    min_len = min(len(ref), len(query))
    mutations: List[PointMutation] = []

    for i in range(min_len):
        r_base = ref[i]
        q_base = query[i]

        if r_base != q_base:
            pos_1based = i + 1
            mut_type = classify_nucleotide_substitution(r_base, q_base)

            # Extraer codón y calcular efecto en la proteína
            ref_codon, codon_pos, offset = get_codon_at_position(ref, i, frame)
            query_codon, _, _ = get_codon_at_position(query, i, frame)

            if ref_codon and query_codon and len(ref_codon) == 3 and len(query_codon) == 3:
                ref_aa = translate_codon(ref_codon)
                alt_aa = translate_codon(query_codon)
                effect = classify_protein_effect(ref_aa, alt_aa)
                desc = (
                    f"Posición {pos_1based}: {r_base} > {q_base} ({mut_type}). "
                    f"Codón {codon_pos}: {ref_codon}({ref_aa}) > {query_codon}({alt_aa}) -> {effect}"
                )
            else:
                ref_aa = None
                alt_aa = None
                effect = "No aplicable (fuera de marco o codón incompleto)"
                desc = f"Posición {pos_1based}: {r_base} > {q_base} ({mut_type})"

            mutations.append(PointMutation(
                position=pos_1based,
                ref_base=r_base,
                alt_base=q_base,
                type=mut_type,
                codon_pos=codon_pos,
                ref_codon=ref_codon,
                alt_codon=query_codon,
                ref_aa=ref_aa,
                alt_aa=alt_aa,
                effect=effect,
                description=desc
            ))

    return mutations
