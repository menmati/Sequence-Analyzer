"""
Módulo para dar formato y generar reportes en texto, consola y JSON.
"""

import json
from typing import Any, Dict, List, Optional
from .codons import CodonMatch, ORF
from .mutations import PointMutation


def format_sequence_display(seq: str, line_width: int = 60, indent: str = "  ") -> str:
    """Divide una secuencia larga en bloques legibles de ancho fijo."""
    lines = []
    for i in range(0, len(seq), line_width):
        chunk = seq[i:i + line_width]
        lines.append(f"{indent}{i+1:5d}  {chunk}")
    return "\n".join(lines)


def format_report_text(analysis: Dict[str, Any]) -> str:
    """Genera un informe completo en texto claro y estructurado."""
    lines: List[str] = []
    sep = "=" * 70
    subsep = "-" * 70

    lines.append(sep)
    lines.append(f" INFORME DE ANÁLISIS DE SECUENCIA: {analysis['id']}")
    if analysis.get('description'):
        lines.append(f" Descripción: {analysis['description']}")
    lines.append(sep)

    # 1. Información General y Contenido GC
    gc = analysis["gc"]
    lines.append("\n[1] PROPIEDADES DE LA SECUENCIA Y CONTENIDO GC")
    lines.append(subsep)
    lines.append(f"  Longitud total:       {analysis['length']} nucleótidos")
    lines.append(f"  Contenido GC:         {gc['gc_percent']}%")
    lines.append(f"  Bases estándar:       A: {gc['counts']['A']} ({gc['frequencies']['A']}%) | "
                 f"T: {gc['counts']['T']} ({gc['frequencies']['T']}%) | "
                 f"C: {gc['counts']['C']} ({gc['frequencies']['C']}%) | "
                 f"G: {gc['counts']['G']} ({gc['frequencies']['G']}%)")
    if gc['counts'].get('otros', 0) > 0:
        lines.append(f"  Bases no estándar/N:  {gc['counts']['otros']}")

    # 2. Codones de Inicio y Parada
    starts: List[CodonMatch] = analysis.get("start_codons", [])
    stops: List[CodonMatch] = analysis.get("stop_codons", [])

    lines.append("\n[2] CODONES DE INICIO Y PARADA")
    lines.append(subsep)
    lines.append(f"  Codones de inicio encontrados (ATG): {len(starts)}")
    if starts:
        start_desc = [f"pos {m.position} (marco +{m.frame})" for m in starts[:8]]
        if len(starts) > 8:
            start_desc.append(f"... y {len(starts) - 8} más")
        lines.append(f"    -> {', '.join(start_desc)}")

    lines.append(f"  Codones de parada encontrados (TAA, TAG, TGA): {len(stops)}")
    if stops:
        stop_desc = [f"{m.codon} en pos {m.position} (marco +{m.frame})" for m in stops[:8]]
        if len(stops) > 8:
            stop_desc.append(f"... y {len(stops) - 8} más")
        lines.append(f"    -> {', '.join(stop_desc)}")

    # ORFs detectados
    orfs: List[ORF] = analysis.get("orfs", [])
    if orfs:
        lines.append(f"\n  Marcos Abiertos de Lectura (ORFs) identificados: {len(orfs)}")
        for idx, orf in enumerate(orfs[:5], 1):
            lines.append(f"    ORF #{idx}: Marco +{orf.frame} | Posición {orf.start_pos}..{orf.stop_pos} "
                         f"| Longitud: {orf.length_nt} nt ({orf.length_aa} aa)")
        if len(orfs) > 5:
            lines.append(f"    ... y {len(orfs) - 5} ORFs adicionales detectados.")

    # 3. Traducción a Proteína
    translations: Dict[str, str] = analysis.get("translations", {})
    lines.append("\n[3] TRADUCCIÓN A PROTEÍNA")
    lines.append(subsep)
    for frame_name, prot in translations.items():
        lines.append(f"  Marco de lectura {frame_name} ({len(prot)} aminoácidos):")
        if len(prot) <= 120:
            lines.append(f"    {prot}")
        else:
            lines.append(format_sequence_display(prot, line_width=60, indent="    "))
        lines.append("")

    # 4. Mutaciones puntuales
    mutations: Optional[List[PointMutation]] = analysis.get("mutations")
    if mutations is not None:
        ref_id = analysis.get("reference_id", "Referencia")
        lines.append(f"\n[4] ANÁLISIS DE MUTACIONES PUNTUALES (frente a '{ref_id}')")
        lines.append(subsep)
        if not mutations:
            lines.append("  No se detectaron mutaciones puntuales. La secuencia es idéntica a la referencia.")
        else:
            lines.append(f"  Total de mutaciones puntuales detectadas: {len(mutations)}")

            # Estadísticas
            transiciones = sum(1 for m in mutations if m.type == "Transición")
            transversiones = sum(1 for m in mutations if m.type == "Transversión")
            ti_tv = f"{transiciones / transversiones:.2f}" if transversiones > 0 else "N/A"
            sinonimas = sum(1 for m in mutations if m.effect and "Sinónima" in m.effect)
            no_sinonimas = sum(1 for m in mutations if m.effect and "No sinónima" in m.effect)
            sin_sentido = sum(1 for m in mutations if m.effect and "Sin sentido" in m.effect)

            lines.append(f"  - Transiciones (Ti): {transiciones} | Transversiones (Tv): {transversiones} | Ratio Ti/Tv: {ti_tv}")
            lines.append(f"  - Efecto funcional: {sinonimas} Sinónimas | {no_sinonimas} Missense | {sin_sentido} Nonsense")
            lines.append("\n  Detalle de sustituciones:")
            lines.append(f"  {'POS':<6} {'REF':<4} {'ALT':<4} {'TIPO':<14} {'CAMBIO CODÓN / AA':<22} {'EFECTO'}")
            lines.append(f"  {'-'*5:<6} {'-'*3:<4} {'-'*3:<4} {'-'*13:<14} {'-'*20:<22} {'-'*20}")

            for m in mutations:
                codon_str = f"{m.ref_codon}({m.ref_aa})>{m.alt_codon}({m.alt_aa})" if m.ref_codon else "-"
                lines.append(f"  {m.position:<6} {m.ref_base:<4} {m.alt_base:<4} {m.type:<14} {codon_str:<22} {m.effect or '-'}")

    lines.append("\n" + sep + "\n")
    return "\n".join(lines)


def format_report_json(analysis_list: List[Dict[str, Any]]) -> str:
    """Convierte los análisis a formato JSON estructurado."""
    serializable = []
    for item in analysis_list:
        data = dict(item)
        # Convertir dataclasses
        if "start_codons" in data:
            data["start_codons"] = [m.__dict__ for m in data["start_codons"]]
        if "stop_codons" in data:
            data["stop_codons"] = [m.__dict__ for m in data["stop_codons"]]
        if "orfs" in data:
            data["orfs"] = [m.__dict__ for m in data["orfs"]]
        if "mutations" in data and data["mutations"] is not None:
            data["mutations"] = [m.__dict__ for m in data["mutations"]]
        serializable.append(data)

    return json.dumps(serializable, indent=2, ensure_ascii=False)
