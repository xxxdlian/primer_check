# Primer Alignment and Quality Check

This script aligns primer sequences to reference sequences with multiple-sequence alignment, extracts the primer binding region with 10 bp flanking sequence on both sides, and exports both FASTA and JPG visualization files for each primer and probe. The primer quality is calculated by `primer3-py` and results can be found in primer_quality.tsv
current version 2


## Input Files

Place these files in the same directory as `primer_check.py`:

- `reference.fasta`: reference FASTA alignment.
- `primers.fasta`: primer sequences in FASTA format.

The reference sequences must already be aligned.

## Primer Naming Rules

Primer names must contain `_`.

The part before `_` is the assay name. Primers with the same assay name are treated as one primer set.

The part after `_` must be one of:

- `forward`
- `reverse`
- `probe`

Examples:

```text
AS_forward
AS_reverse
AS_probe
BA_forward
BA_reverse
BA_probe
```

Each assay must contain one `forward` primer and one `reverse` primer. A `probe` is optional.

## Reverse Primer Handling

Reverse primers use the `reverse` suffix. The script reverse-complements that primer sequence before alignment and quality calculation.

Examples:

- `AS_reverse` is reverse-complemented.
- `AS_forward` is used as-is.

## Installation

Create or activate your conda environment, then install the required packages:

```bash
conda activate your_env_name
pip install pillow primer3-py
```

`pillow` is required for JPG output.

`primer3-py` is needed for quality check

## Run

From this project directory:

```bash
python3 primer_check.py
```

## Output

All output files are written to:

```text
primer_msa_fastas/
```

The script creates:

- One local FASTA alignment per primer:

```text
primer_msa_fastas/RSU_forward.fasta
primer_msa_fastas/RSU_reverse.fasta
...
```

- One JPG visualization per primer:

```text
primer_msa_fastas/jpg/RSU_forward.jpg
primer_msa_fastas/jpg/RSU_reverse.jpg
...
```

- One probe map JPG per assay when a probe is present:

```text
primer_msa_fastas/probe_maps/AS_probe_map.jpg
primer_msa_fastas/probe_maps/BA_probe_map.jpg
```

- One primer quality summary table:

```text
primer_msa_fastas/primer_quality.tsv
```


## FASTA Output

Each primer FASTA contains:

- all reference sequences in the local primer region
- the primer sequence aligned to that same region
- 10 bp of flanking sequence on each side of the primer binding site

The primer row keeps the full primer sequence. Positions outside the primer are padded with `-`.

## JPG Output

Each JPG shows the local alignment region.

Display rules:

- The primer sequence is shown first and is displayed in full.
- Reference sequences are compared against the primer sequence.
- Reference positions matching the primer are shown as `-`.
- Reference positions different from the primer are shown as actual bases.

## Primer Quality TSV Columns

`primer_quality.tsv` includes:

- primer name
- sequence used for analysis
- whether the sequence was original or reverse-complemented
- length
- GC percentage
- Primer3 Tm
- hairpin Tm and delta G
- self-dimer Tm and delta G
- forward/reverse hetero-dimer Tm and delta G
- PCR product length
- probe name and probe position inside the PCR product, when present
- probe map JPG path, when present
- best matching reference
- primer binding coordinates
- extracted local alignment coordinates
- FASTA output path
- JPG output path


## Notes

Coordinates in the output table are 1-based.

The script uses ungapped reference coordinates to find primer positions, then maps those positions back to the original gapped alignment columns for FASTA and JPG output.
