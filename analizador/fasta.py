"""
Módulo para el procesamiento y lectura de archivos en formato FASTA.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import List, Union


@dataclass
class FastaRecord:
    """Representa una secuencia individual dentro de un archivo FASTA."""
    id: str
    description: str
    sequence: str

    @property
    def length(self) -> int:
        return len(self.sequence)

    def __repr__(self) -> str:
        return f"<FastaRecord id='{self.id}' len={self.length}>"


def parse_fasta(source: Union[str, Path]) -> List[FastaRecord]:
    """
    Lee un archivo FASTA (o contenido de texto) y devuelve una lista de FastaRecord.

    Args:
        source: Ruta al archivo FASTA o cadena de texto con formato FASTA.

    Returns:
        Lista de instancias de FastaRecord.

    Raises:
        FileNotFoundError: Si el archivo especificado no existe.
        ValueError: Si el formato FASTA es inválido.
    """
    path = Path(source)
    if path.is_file():
        with open(path, "r", encoding="utf-8") as f:
            lines = f.readlines()
    elif isinstance(source, str) and (">" in source):
        lines = source.splitlines(keepends=True)
    else:
        raise FileNotFoundError(f"No se encontró el archivo: {source}")

    records: List[FastaRecord] = []
    current_id = ""
    current_desc = ""
    current_seq_parts: List[str] = []

    for line_num, line in enumerate(lines, start=1):
        line = line.strip()
        if not line:
            continue

        if line.startswith(">"):
            if current_id or current_seq_parts:
                full_seq = "".join(current_seq_parts).upper()
                records.append(FastaRecord(id=current_id, description=current_desc, sequence=full_seq))
                current_seq_parts = []

            header = line[1:].strip()
            parts = header.split(maxsplit=1)
            current_id = parts[0] if parts else f"seq_{len(records) + 1}"
            current_desc = parts[1] if len(parts) > 1 else ""
        else:
            # Línea de secuencia: limpiar espacios y validar caracteres básicos
            cleaned = "".join(c for c in line if not c.isspace()).upper()
            current_seq_parts.append(cleaned)

    if current_id or current_seq_parts:
        full_seq = "".join(current_seq_parts).upper()
        records.append(FastaRecord(id=current_id, description=current_desc, sequence=full_seq))

    if not records:
        raise ValueError("El archivo no contiene secuencias en formato FASTA válido.")

    return records


COMPLEMENT_MAP = str.maketrans({
    "A": "T", "T": "A", "U": "A",
    "C": "G", "G": "C",
    "R": "Y", "Y": "R",
    "S": "S", "W": "W",
    "K": "M", "M": "K",
    "B": "V", "V": "B",
    "D": "H", "H": "D",
    "N": "N"
})


def reverse_complement(sequence: str) -> str:
    """
    Calcula la hebra reversa complementaria de una secuencia de ADN/ARN.
    """
    clean_seq = sequence.upper()
    return clean_seq.translate(COMPLEMENT_MAP)[::-1]
