#!/usr/bin/env python3
from pathlib import Path
import re

from PIL import Image, ImageDraw, ImageFont

try:
    import primer3
except ImportError:
    primer3 = None


REFERENCE_FASTA = Path("reference.fasta")
PRIMER_FASTA = Path("primers.fasta")
OUT_DIR = Path("primer_msa_fastas")
FLANK_BP = 10
JPEG_DIR = OUT_DIR / "jpg"
PROBE_MAP_DIR = OUT_DIR / "probe_maps"
QUALITY_TSV = OUT_DIR / "primer_quality.tsv"
VALID_PRIMER_ROLES = {"forward", "reverse", "probe"}


def read_fasta(path):
    records = []
    name = None
    desc = None
    seq = []
    with path.open() as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append((name, desc, "".join(seq).upper()))
                desc = line[1:]
                name = desc.split()[0]
                seq = []
            else:
                seq.append(line)
    if name is not None:
        records.append((name, desc, "".join(seq).upper()))
    return records


def safe_name(name):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", name)


def parse_primer_name(name):
    if "_" not in name:
        raise ValueError(
            f"Primer name '{name}' is invalid. Expected format: assay_forward, assay_reverse, or assay_probe."
        )

    assay, role = name.split("_", 1)
    role = role.lower()
    if not assay or role not in VALID_PRIMER_ROLES:
        raise ValueError(
            f"Primer name '{name}' is invalid. The suffix after '_' must be forward, reverse, or probe."
        )

    return assay, role


def reverse_complement(seq):
    table = str.maketrans(
        "ACGTRYKMSWBDHVNacgtrykmswbdhvn",
        "TGCAYRMKSWVHDBNtgcayrmkswvhdbn",
    )
    return seq.translate(table)[::-1].upper()


def gc_percent(seq):
    if not seq:
        return 0
    return 100 * sum(base in "GC" for base in seq.upper()) / len(seq)


def primer3_quality(seq):
    if primer3 is None:
        return {
            "tm": "NA",
            "gc_percent": f"{gc_percent(seq):.2f}",
            "hairpin_tm": "NA",
            "hairpin_dg": "NA",
            "hairpin_structure_found": "NA",
            "self_dimer_tm": "NA",
            "self_dimer_dg": "NA",
            "self_dimer_structure_found": "NA",
        }

    hairpin = primer3.calc_hairpin(seq)
    self_dimer = primer3.calc_homodimer(seq)
    return {
        "tm": f"{primer3.calc_tm(seq):.2f}",
        "gc_percent": f"{gc_percent(seq):.2f}",
        "hairpin_tm": f"{hairpin.tm:.2f}",
        "hairpin_dg": f"{hairpin.dg / 1000:.2f}",
        "hairpin_structure_found": str(hairpin.structure_found),
        "self_dimer_tm": f"{self_dimer.tm:.2f}",
        "self_dimer_dg": f"{self_dimer.dg / 1000:.2f}",
        "self_dimer_structure_found": str(self_dimer.structure_found),
    }


def heterodimer_quality(forward_seq, reverse_seq):
    if primer3 is None:
        return {
            "pair_heterodimer_tm": "NA",
            "pair_heterodimer_dg": "NA",
            "pair_heterodimer_structure_found": "NA",
        }

    heterodimer = primer3.calc_heterodimer(forward_seq, reverse_seq)
    return {
        "pair_heterodimer_tm": f"{heterodimer.tm:.2f}",
        "pair_heterodimer_dg": f"{heterodimer.dg / 1000:.2f}",
        "pair_heterodimer_structure_found": str(heterodimer.structure_found),
    }


def best_primer_position(primer_seq, references):
    best = None
    for ref_name, ref_desc, ref_seq in references:
        ungapped_cols = [i for i, base in enumerate(ref_seq) if base != "-"]
        ungapped_seq = "".join(ref_seq[i] for i in ungapped_cols)
        if len(ungapped_seq) < len(primer_seq):
            continue

        for start in range(len(ungapped_seq) - len(primer_seq) + 1):
            window = ungapped_seq[start : start + len(primer_seq)]
            current_matches = sum(p == r for p, r in zip(primer_seq, window))
            if best is None or current_matches > best["matches"]:
                best = {
                    "matches": current_matches,
                    "ref_name": ref_name,
                    "ref_desc": ref_desc,
                    "ungapped_start": start,
                    "ungapped_end": start + len(primer_seq) - 1,
                    "cols": ungapped_cols[start : start + len(primer_seq)],
                }
    return best


def region_with_flank(best, references, flank_bp):
    for ref_name, _, ref_seq in references:
        if ref_name != best["ref_name"]:
            continue

        ungapped_cols = [i for i, base in enumerate(ref_seq) if base != "-"]
        region_start = max(0, best["ungapped_start"] - flank_bp)
        region_end = min(len(ungapped_cols) - 1, best["ungapped_end"] + flank_bp)
        return ungapped_cols[region_start], ungapped_cols[region_end]

    raise ValueError(f"Best reference not found: {best['ref_name']}")


def wrap(seq, width=80):
    return "\n".join(seq[i : i + width] for i in range(0, len(seq), width))


def short_label(desc, max_len=34):
    label = desc.split()[0]
    return label[:max_len]


def render_alignment_jpg(records, primer_name, out_path):
    font = ImageFont.load_default()
    left_margin = 250
    top_margin = 32
    row_height = 22
    char_width = 10
    char_height = 16
    label_color = (35, 35, 35)
    bg = (255, 255, 255)
    colors = {
        "A": (131, 198, 136),
        "C": (121, 169, 218),
        "G": (246, 198, 103),
        "T": (230, 128, 120),
        "U": (230, 128, 120),
        "-": (235, 235, 235),
    }

    aln_len = max(len(seq) for _, seq in records)
    width = left_margin + aln_len * char_width + 24
    height = top_margin + len(records) * row_height + 24
    image = Image.new("RGB", (width, height), bg)
    draw = ImageDraw.Draw(image)
    draw.text((12, 8), primer_name, fill=label_color, font=font)

    reference_seq = records[0][1] if records else ""

    for row_idx, (label, seq) in enumerate(records):
        y = top_margin + row_idx * row_height
        draw.text((12, y + 3), short_label(label), fill=label_color, font=font)
        for col_idx, base in enumerate(seq):
            display_base = base
            if row_idx > 0 and col_idx < len(reference_seq):
                if base.upper() == reference_seq[col_idx].upper():
                    display_base = "-"

            x = left_margin + col_idx * char_width
            color = colors.get(display_base.upper(), colors.get(base.upper(), (210, 210, 210)))
            draw.rectangle(
                [x, y, x + char_width - 1, y + char_height],
                fill=color,
                outline=(255, 255, 255),
            )
            draw.text((x + 2, y + 3), display_base, fill=(25, 25, 25), font=font)

    image.save(out_path, quality=95)


def render_probe_map_jpg(assay, forward_result, reverse_result, probe_result, out_path):
    font = ImageFont.load_default()
    margin = 48
    width = 900
    height = 230
    line_y = 104
    product_start = min(
        forward_result["best"]["ungapped_start"],
        reverse_result["best"]["ungapped_start"],
    )
    product_end = max(
        forward_result["best"]["ungapped_end"],
        reverse_result["best"]["ungapped_end"],
    )
    product_len = product_end - product_start + 1

    def x_from_pos(pos):
        if product_len <= 1:
            return margin
        return margin + int((pos - product_start) / (product_len - 1) * (width - 2 * margin))

    image = Image.new("RGB", (width, height), (255, 255, 255))
    draw = ImageDraw.Draw(image)
    draw.text((margin, 18), f"{assay} PCR product", fill=(30, 30, 30), font=font)
    draw.text(
        (margin, 36),
        f"product length: {product_len} bp",
        fill=(70, 70, 70),
        font=font,
    )
    draw.line((margin, line_y, width - margin, line_y), fill=(55, 55, 55), width=4)

    blocks = [
        ("forward", forward_result, (96, 167, 102), line_y - 34),
        ("probe", probe_result, (238, 187, 83), line_y - 8),
        ("reverse", reverse_result, (213, 103, 99), line_y + 18),
    ]
    for label, result, color, y in blocks:
        start = result["best"]["ungapped_start"]
        end = result["best"]["ungapped_end"]
        x1 = x_from_pos(start)
        x2 = max(x1 + 6, x_from_pos(end))
        draw.rectangle((x1, y, x2, y + 18), fill=color, outline=(35, 35, 35))
        draw.text((x1, y - 14), label, fill=(30, 30, 30), font=font)
        draw.text(
            (x1, y + 22),
            f"{start + 1}-{end + 1}",
            fill=(70, 70, 70),
            font=font,
        )

    draw.text((margin, height - 34), f"start {product_start + 1}", fill=(70, 70, 70), font=font)
    draw.text(
        (width - margin - 70,
        height - 34),
        f"end {product_end + 1}",
        fill=(70, 70, 70),
        font=font,
    )
    image.save(out_path, quality=95)


def main():
    references = read_fasta(REFERENCE_FASTA)
    primers = read_fasta(PRIMER_FASTA)
    if not references:
        raise SystemExit(f"No records found in {REFERENCE_FASTA}")
    if not primers:
        raise SystemExit(f"No records found in {PRIMER_FASTA}")

    aln_len = len(references[0][2])
    for ref_name, _, ref_seq in references:
        if len(ref_seq) != aln_len:
            raise SystemExit(f"Reference alignment lengths differ at {ref_name}")

    OUT_DIR.mkdir(exist_ok=True)
    JPEG_DIR.mkdir(exist_ok=True)
    PROBE_MAP_DIR.mkdir(exist_ok=True)

    grouped_primers = {}
    for primer_name, _, _ in primers:
        assay, role = parse_primer_name(primer_name)
        grouped_primers.setdefault(assay, {})
        if role in grouped_primers[assay]:
            raise SystemExit(f"Duplicate {role} primer found for assay '{assay}'.")
        grouped_primers[assay][role] = primer_name

    for assay, assay_primers in grouped_primers.items():
        if "forward" not in assay_primers or "reverse" not in assay_primers:
            raise SystemExit(f"Assay '{assay}' must contain both forward and reverse primers.")

    quality_header = [
        "assay",
        "role",
        "primer",
        "sequence_used",
        "sequence_analyzed",
        "length",
        "gc_percent",
        "tm",
        "hairpin_tm",
        "hairpin_dg_kcal_mol",
        "hairpin_structure_found",
        "self_dimer_tm",
        "self_dimer_dg_kcal_mol",
        "self_dimer_structure_found",
        "pair_forward",
        "pair_reverse",
        "pair_heterodimer_tm",
        "pair_heterodimer_dg_kcal_mol",
        "pair_heterodimer_structure_found",
        "pcr_product_length_bp",
        "probe_name",
        "probe_product_start",
        "probe_product_end",
        "probe_map_file",
        "best_reference",
        "ungapped_start",
        "ungapped_end",
        "alignment_start",
        "alignment_end",
        "region_alignment_start",
        "region_alignment_end",
        "flank_bp",
        "fasta_file",
        "jpg_file",
    ]
    quality_lines = ["\t".join(quality_header)]
    results = {}

    for primer_name, primer_desc, primer_seq in primers:
        assay, role = parse_primer_name(primer_name)
        primer_seq = primer_seq.replace("-", "")
        sequence_used = "original"
        if role == "reverse":
            primer_seq = reverse_complement(primer_seq)
            sequence_used = "reverse_complement"
        quality = primer3_quality(primer_seq)

        best = best_primer_position(primer_seq, references)
        if best is None:
            print(f"Skip {primer_name}: no valid position found")
            continue

        primer_aln = ["-"] * aln_len
        for base, col in zip(primer_seq, best["cols"]):
            primer_aln[col] = base
        primer_aln = "".join(primer_aln)
        region_start_col, region_end_col = region_with_flank(best, references, FLANK_BP)

        out_path = OUT_DIR / f"{safe_name(primer_name)}.fasta"
        jpg_path = JPEG_DIR / f"{safe_name(primer_name)}.jpg"
        local_records = []
        primer_local_seq = primer_aln[region_start_col : region_end_col + 1]
        with out_path.open("w") as out:
            for _, ref_desc, ref_seq in references:
                local_seq = ref_seq[region_start_col : region_end_col + 1]
                local_records.append((ref_desc, local_seq))
                out.write(
                    f">{ref_desc} region_alignment_cols={region_start_col + 1}-{region_end_col + 1}\n"
                    f"{wrap(local_seq)}\n"
                )
            local_records.append((primer_desc, primer_local_seq))
            out.write(
                f">{primer_desc} aligned_to={best['ref_name']} "
                f"sequence_used={sequence_used} "
                f"ungapped={best['ungapped_start'] + 1}-{best['ungapped_end'] + 1} "
                f"alignment_cols={best['cols'][0] + 1}-{best['cols'][-1] + 1} "
                f"region_alignment_cols={region_start_col + 1}-{region_end_col + 1} "
                f"flank_bp={FLANK_BP}\n"
            )
            out.write(f"{wrap(primer_local_seq)}\n")

        jpg_records = [(primer_desc, primer_local_seq)]
        for _, ref_desc, ref_seq in references:
            local_seq = ref_seq[region_start_col : region_end_col + 1]
            jpg_records.append((ref_desc, local_seq))
        render_alignment_jpg(jpg_records, primer_name, jpg_path)

        results[primer_name] = {
            "assay": assay,
            "role": role,
            "primer_name": primer_name,
            "sequence_used": sequence_used,
            "primer_seq": primer_seq,
            "quality": quality,
            "best": best,
            "region_start_col": region_start_col,
            "region_end_col": region_end_col,
            "out_path": out_path,
            "jpg_path": jpg_path,
        }

    assay_metrics = {}
    for assay, assay_primers in grouped_primers.items():
        forward = results[assay_primers["forward"]]
        reverse = results[assay_primers["reverse"]]
        product_start = min(
            forward["best"]["ungapped_start"],
            reverse["best"]["ungapped_start"],
        )
        product_end = max(
            forward["best"]["ungapped_end"],
            reverse["best"]["ungapped_end"],
        )
        product_length = product_end - product_start + 1
        pair_quality = heterodimer_quality(forward["primer_seq"], reverse["primer_seq"])

        probe_name = assay_primers.get("probe", "NA")
        probe_start = "NA"
        probe_end = "NA"
        probe_map_file = "NA"
        if probe_name != "NA":
            probe = results[probe_name]
            probe_start = str(probe["best"]["ungapped_start"] - product_start + 1)
            probe_end = str(probe["best"]["ungapped_end"] - product_start + 1)
            probe_map_path = PROBE_MAP_DIR / f"{safe_name(assay)}_probe_map.jpg"
            render_probe_map_jpg(assay, forward, reverse, probe, probe_map_path)
            probe_map_file = str(probe_map_path)

        assay_metrics[assay] = {
            "pair_forward": forward["primer_name"],
            "pair_reverse": reverse["primer_name"],
            "pair_quality": pair_quality,
            "pcr_product_length": str(product_length),
            "probe_name": probe_name,
            "probe_product_start": probe_start,
            "probe_product_end": probe_end,
            "probe_map_file": probe_map_file,
        }

    for primer_name, _, _ in primers:
        result = results[primer_name]
        quality = result["quality"]
        best = result["best"]
        metrics = assay_metrics[result["assay"]]
        pair_quality = metrics["pair_quality"]
        quality_lines.append(
            "\t".join(
                [
                    result["assay"],
                    result["role"],
                    primer_name,
                    result["sequence_used"],
                    result["primer_seq"],
                    str(len(result["primer_seq"])),
                    quality["gc_percent"],
                    quality["tm"],
                    quality["hairpin_tm"],
                    quality["hairpin_dg"],
                    quality["hairpin_structure_found"],
                    quality["self_dimer_tm"],
                    quality["self_dimer_dg"],
                    quality["self_dimer_structure_found"],
                    metrics["pair_forward"],
                    metrics["pair_reverse"],
                    pair_quality["pair_heterodimer_tm"],
                    pair_quality["pair_heterodimer_dg"],
                    pair_quality["pair_heterodimer_structure_found"],
                    metrics["pcr_product_length"],
                    metrics["probe_name"],
                    metrics["probe_product_start"],
                    metrics["probe_product_end"],
                    metrics["probe_map_file"],
                    best["ref_name"],
                    str(best["ungapped_start"] + 1),
                    str(best["ungapped_end"] + 1),
                    str(best["cols"][0] + 1),
                    str(best["cols"][-1] + 1),
                    str(result["region_start_col"] + 1),
                    str(result["region_end_col"] + 1),
                    str(FLANK_BP),
                    str(result["out_path"]),
                    str(result["jpg_path"]),
                ]
            )
        )

    QUALITY_TSV.write_text("\n".join(quality_lines) + "\n")

    print(f"Wrote {len(primers)} primer MSA FASTA files to {OUT_DIR}")
    print(f"Wrote alignment JPG files to {JPEG_DIR}")
    print(f"Wrote probe map JPG files to {PROBE_MAP_DIR}")
    print(f"Primer quality TSV: {QUALITY_TSV}")
    if primer3 is None:
        print("Warning: primer3-py is not installed; primer3 quality fields were written as NA.")


if __name__ == "__main__":
    main()
