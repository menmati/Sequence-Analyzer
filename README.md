# Analizador de Secuencias FASTA

Herramienta de línea de comandos (CLI) desarrollada en Python para el análisis bioinformático de secuencias de ADN y ARN en formato FASTA.

## Características

1. **Lectura y Validación de Archivos FASTA**:
   - Admite secuencias simples y multilínea, múltiples registros por archivo y caracteres estándar/ambiguos.
2. **Cálculo del Contenido de GC**:
   - Porcentaje global de GC.
   - Frecuencias y conteo individual de bases (A, C, G, T y ambiguas).
3. **Identificación de Codones de Inicio y Parada**:
   - Detección de codones de inicio (`ATG`) y codones de parada (`TAA`, `TAG`, `TGA`).
   - Ubicación de posiciones exactas (en base 1) y marcos de lectura asociados.
   - Detección de Marcos Abiertos de Lectura (**ORFs**).
4. **Traducción a Secuencias Proteicas**:
   - Código genético estándar universal.
   - Traducción en marcos directos (`+1`, `+2`, `+3`) y opcionalmente reversos complementarios (`-1`, `-2`, `-3`).
5. **Identificación y Clasificación de Mutaciones Puntuales**:
   - Detección de sustituciones nucleotídicas frente a una secuencia de referencia.
   - Clasificación en **Transición** (purina↔purina o pirimidina↔pirimidina) o **Transversión** (purina↔pirimidina).
   - Determinación del impacto funcional en la proteína:
     - **Sinónima (Silenciosa)**: Mismo aminoácido codificado.
     - **No sinónima (Missense / Cambio de sentido)**: Cambio a otro aminoácido.
     - **Sin sentido (Nonsense)**: Generación de codón de parada prematuro.
     - **Pérdida de parada (Readthrough)**.
   - Cálculo del ratio Ti/Tv.
6. **Sin dependencias externas**:
   - Funciona directamente con la biblioteca estándar de Python (3.8+).
   - Compatible con exportación en formato texto legible o **JSON**.

---

## Estructura del Proyecto

```
analizador_secuencias/
├── analizador/
│   ├── __init__.py       # Exportaciones principales del paquete
│   ├── fasta.py          # Lector y utilidades de formato FASTA
│   ├── codons.py         # Cálculo de GC, detección de codones y ORFs
│   ├── translation.py    # Código genético y traducción a proteínas
│   ├── mutations.py      # Detección y clasificación de mutaciones puntuales
│   └── report.py         # Generador de reportes en consola y JSON
├── ejemplos/
│   ├── ejemplo_multiples.fasta  # Gen de beta-globina (wt y variantes)
│   ├── referencia.fasta         # Secuencia de referencia de ejemplo
│   └── mutante.fasta            # Variante con mutaciones de prueba
├── tests/
│   └── test_analizador.py       # Pruebas unitarias completas
├── main.py               # Punto de entrada de la línea de comandos
└── README.md             # Documentación de uso
```

---

## Uso de la Herramienta

### 1. Análisis básico de un archivo FASTA

```bash
python main.py ejemplos/ejemplo_multiples.fasta
```

### 2. Comparación de mutaciones frente a un archivo de referencia externo

```bash
python main.py ejemplos/mutante.fasta --referencia ejemplos/referencia.fasta
```

### 3. Seleccionar marco de lectura específico

Traducir únicamente en el marco de lectura 1:

```bash
python main.py ejemplos/referencia.fasta --marco 1
```

Opciones disponibles para `--marco`: `1`, `2`, `3` o `todos` (por defecto).

### 4. Traducir ambas hebras (directa y reversa complementaria)

```bash
python main.py ejemplos/referencia.fasta --ambas-hebras
```

### 5. Exportar el informe a formato JSON

```bash
python main.py ejemplos/mutante.fasta --referencia ejemplos/referencia.fasta --formato json --salida resultado.json
```

---

## Opciones de la Línea de Comandos

| Parámetro | Argumento | Descripción |
|---|---|---|
| `archivo_fasta` | Posicional | Ruta al archivo FASTA a analizar. |
| `-r`, `--referencia` | `<ruta o ID>` | Archivo FASTA o ID de secuencia de referencia para mutaciones. |
| `-m`, `--marco` | `1`, `2`, `3`, `todos` | Marco de lectura para la traducción (por defecto: `todos`). |
| `--min-orf` | `<entero>` | Longitud mínima de nucleótidos para considerar un ORF (defecto: 30 nt). |
| `--ambas-hebras` | Flag | Traduce marcos `+1, +2, +3` y `-1, -2, -3` de la hebra complementaria. |
| `-f`, `--formato` | `texto`, `json` | Formato del informe de salida (por defecto: `texto`). |
| `-o`, `--salida` | `<archivo>` | Guarda el reporte en un archivo en vez de la consola. |

---

## Ejecución de Pruebas Unitarias

Para ejecutar la suite de pruebas automatizadas:

```bash
python -m unittest tests/test_analizador.py
```
