#!/usr/bin/env python3
"""
Herramienta de Línea de Comandos (CLI) para Análisis de Secuencias FASTA.

Funcionalidades:
- Lectura y validación de archivos FASTA
- Cálculo del contenido de GC y proporciones de nucleótidos
- Identificación de codones de inicio y parada
- Detección de Marcos Abiertos de Lectura (ORFs)
- Traducción a secuencias proteicas en diferentes marcos de lectura
- Identificación y caracterización de mutaciones puntuales (SNPs, transiciones,
  transversiones, sinónimas, no sinónimas y codones de parada prematuros)
"""

import argparse
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

# Asegurar codificación UTF-8 en consola para caracteres en español
if sys.stdout.encoding and sys.stdout.encoding.lower() != "utf-8":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
        sys.stderr.reconfigure(encoding="utf-8")
    except Exception:
        pass

from analizador.fasta import parse_fasta, FastaRecord
from analizador.codons import calculate_gc_content, find_codons, find_orfs
from analizador.translation import translate_dna, translate_all_frames
from analizador.mutations import identify_point_mutations, PointMutation
from analizador.report import format_report_text, format_report_json


def analyze_sequence(
    record: FastaRecord,
    frame_choice: str = "todos",
    min_orf_len: int = 30,
    both_strands: bool = False,
    ref_record: Optional[FastaRecord] = None
) -> Dict[str, Any]:
    """Ejecuta el análisis biológico completo para una secuencia FASTA."""
    seq = record.sequence
    gc_data = calculate_gc_content(seq)
    starts, stops = find_codons(seq)
    orfs = find_orfs(seq, min_length_nt=min_orf_len)

    # Traducción
    if frame_choice == "todos":
        translations = translate_all_frames(seq, include_reverse=both_strands)
    else:
        frame_int = int(frame_choice)
        translations = {f"+{frame_int}": translate_dna(seq, frame=frame_int)}

    # Análisis de mutaciones si se proporciona referencia
    mutations: Optional[List[PointMutation]] = None
    ref_id = None
    if ref_record is not None and ref_record.id != record.id:
        ref_id = ref_record.id
        frame_for_mut = int(frame_choice) if frame_choice in {"1", "2", "3"} else 1
        mutations = identify_point_mutations(ref_record.sequence, seq, frame=frame_for_mut)

    return {
        "id": record.id,
        "description": record.description,
        "length": record.length,
        "gc": gc_data,
        "start_codons": starts,
        "stop_codons": stops,
        "orfs": orfs,
        "translations": translations,
        "mutations": mutations,
        "reference_id": ref_id,
    }


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="analizador-secuencias",
        description="Analizador de secuencias de ADN/ARN en formato FASTA: "
                    "Contenido GC, codones de inicio/parada, traducción proteica y mutaciones puntuales.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Ejemplos de uso:
  # Análisis general de un archivo FASTA
  python main.py secuencia.fasta

  # Comparar secuencias mutantes contra una referencia externa
  python main.py mutantes.fasta --referencia referencia.fasta

  # Traducir sólo el marco de lectura 1 y exportar a JSON
  python main.py gen.fasta --marco 1 --formato json --salida resultado.json

  # Traducir ambas hebras (+1,+2,+3 y -1,-2,-3)
  python main.py secuencia.fasta --ambas-hebras
        """
    )

    parser.add_argument(
        "archivo_fasta",
        type=str,
        help="Ruta al archivo FASTA a analizar."
    )

    parser.add_argument(
        "-r", "--referencia",
        type=str,
        default=None,
        help="Archivo FASTA o ID de secuencia de referencia para detectar mutaciones puntuales. "
             "Si no se especifica y el archivo contiene varias secuencias, la primera se usará como referencia."
    )

    parser.add_argument(
        "-m", "--marco",
        choices=["1", "2", "3", "todos"],
        default="todos",
        help="Marco de lectura para la traducción (1, 2, 3 o 'todos'). Por defecto: 'todos'."
    )

    parser.add_argument(
        "--min-orf",
        type=int,
        default=30,
        help="Longitud mínima en nucleótidos para detectar un ORF (por defecto: 30 nt = 10 aa)."
    )

    parser.add_argument(
        "--ambas-hebras",
        action="store_true",
        help="Traducir también los marcos de lectura de la hebra complementaria reversa (-1, -2, -3)."
    )

    parser.add_argument(
        "-f", "--formato",
        choices=["texto", "json"],
        default="texto",
        help="Formato de salida: 'texto' o 'json' (por defecto: 'texto')."
    )

    parser.add_argument(
        "-o", "--salida",
        type=str,
        default=None,
        help="Ruta al archivo de salida para guardar los resultados (si no se especifica, se muestra en pantalla)."
    )

    return parser.parse_args()


def main():
    args = parse_arguments()

    try:
        query_records = parse_fasta(args.archivo_fasta)
    except Exception as e:
        sys.stderr.write(f"Error al leer el archivo FASTA '{args.archivo_fasta}': {e}\n")
        sys.exit(1)

    ref_record: Optional[FastaRecord] = None

    if args.referencia:
        ref_path = Path(args.referencia)
        if ref_path.is_file():
            try:
                ref_records = parse_fasta(ref_path)
                ref_record = ref_records[0]
            except Exception as e:
                sys.stderr.write(f"Error al leer el archivo de referencia '{args.referencia}': {e}\n")
                sys.exit(1)
        else:
            # Buscar en el mismo archivo por ID
            matching = [r for r in query_records if r.id == args.referencia]
            if matching:
                ref_record = matching[0]
            else:
                sys.stderr.write(f"No se encontró la referencia '{args.referencia}' como archivo ni como ID en el FASTA.\n")
                sys.exit(1)
    elif len(query_records) > 1:
        # Por defecto, si hay varias secuencias en el archivo FASTA, la primera actúa como referencia
        ref_record = query_records[0]

    all_analyses: List[Dict[str, Any]] = []

    for record in query_records:
        analysis = analyze_sequence(
            record=record,
            frame_choice=args.marco,
            min_orf_len=args.min_orf,
            both_strands=args.ambas_hebras,
            ref_record=ref_record
        )
        all_analyses.append(analysis)

    # Generar salida
    if args.formato == "json":
        output_content = format_report_json(all_analyses)
    else:
        text_reports = [format_report_text(a) for a in all_analyses]
        output_content = "\n".join(text_reports)

    if args.salida:
        try:
            with open(args.salida, "w", encoding="utf-8") as f:
                f.write(output_content)
            print(f"Reporte generado exitosamente en: {args.salida}")
        except Exception as e:
            sys.stderr.write(f"Error al guardar el archivo de salida '{args.salida}': {e}\n")
            sys.exit(1)
    else:
        print(output_content)


if __name__ == "__main__":
    main()
