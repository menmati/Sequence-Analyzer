"""
Módulo para análisis de secuencias biológicas en formato FASTA.
"""

from .fasta import FastaRecord, parse_fasta, reverse_complement
from .codons import calculate_gc_content, find_codons, find_orfs, CodonMatch, ORF
from .translation import translate_dna, translate_all_frames, GENETIC_CODE
from .mutations import identify_point_mutations, PointMutation

__all__ = [
    "FastaRecord",
    "parse_fasta",
    "reverse_complement",
    "calculate_gc_content",
    "find_codons",
    "find_orfs",
    "CodonMatch",
    "ORF",
    "translate_dna",
    "translate_all_frames",
    "GENETIC_CODE",
    "identify_point_mutations",
    "PointMutation",
]
