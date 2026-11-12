#!/usr/bin/env python3
"""
prepare_release.py - get the genome-wide splice hotspot output ready for
(1) the GitHub Pages browser and (2) a Zenodo release.

Three steps, run in order:

  1. inventory  Optional. Prefer a hand-written config.json for this repo:
                gene_hotspot_summary.csv + shortlist_interior_hotspots.csv.
                A full inventory walk rescans hundreds of thousands of files
                and can guess the wrong tables.

        python scripts/prepare_release.py inventory --root genome_wide

  2. web        Build the small JSON/CSV files the website loads.

        python scripts/prepare_release.py web --config config.json --out docs/data --workers 8

  3. zenodo     Build the release: summary-table zip, one zip per chromosome
                with per-gene tables and PNGs (SVGs optional), a manifest
                with SHA-256 per file, MD5 per zip, not_scored.tsv and a README.

        python scripts/prepare_release.py zenodo --config config.json --out release

Only the Python standard library is needed (Python 3.8+).  Input files are
only read, never modified.

Web builds use the scored RefSeq transcript (transcript.json + ncbiRefSeq +
*_introns.tsv), not a GENCODE GTF. Zenodo packs only genes present in the
summary table and only the result CSVs listed in config.results_files.
"""
import argparse, csv, fnmatch, glob, gzip, hashlib, json, os, re, sys, zipfile
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import date

csv.field_size_limit(min(sys.maxsize, 2**31 - 1))

INTRON_NAME_RE = re.compile(r"^intron_(\d+)$", re.IGNORECASE)
DEFAULT_RESULTS_FILES = (
    "gene_hotspot_summary.csv",
    "shortlist_interior_hotspots.csv",
    "all_genes_hotspot_intervals.csv",
    "cohort_intron_summary.csv",
)

# Candidate column names used to guess the config; first match wins (case-insensitive).
GUESS = {
    "gene": ["gene", "gene_name", "gene_symbol", "symbol", "hgnc_symbol", "genename"],
    "chrom": ["chrom", "chr", "chromosome", "seqnames", "seqname", "#chrom", "contig"],
    "start": ["start", "hotspot_start", "start_pos", "begin", "chromstart"],
    "end": ["end", "hotspot_end", "end_pos", "stop", "chromend"],
    "intron": ["intron", "intron_number", "intron_id", "intron_idx", "intron_index"],
    "intron_start": ["intron_start"],
    "intron_end": ["intron_end"],
    "strand": ["strand"],
    "length": ["length", "len", "width", "hotspot_length", "size_bp"],
    "q": ["bh_q", "q", "qvalue", "q_value", "padj", "fdr", "bh_qvalue", "p_adj"],
    "n_spliceogenic": ["n_spliceogenic", "spliceogenic_count", "n_splice", "spliceogenic", "n_sig"],
    "max_score": ["max_atlas_score", "max_atlas", "atlas_max", "max_score", "max_spliceai"],
    "pos": ["pos", "position", "start", "bp", "variant_pos"],
    "ref": ["ref", "reference"],
    "alt": ["alt", "alternate"],
    "score": ["atlas_score", "score", "ds_max", "delta_score", "spliceai", "max_delta"],
}


# ---------------------------------------------------------------- helpers
def smart_open(path):
    return gzip.open(path, "rt", newline="") if str(path).endswith(".gz") else open(path, newline="")


def sniff_reader(fh):
    first = fh.readline()
    delim = "\t" if first.count("\t") >= first.count(",") else ","
    header = next(csv.reader([first], delimiter=delim))
    return delim, [h.strip() for h in header], csv.reader(fh, delimiter=delim)


def read_table(path):
    with smart_open(path) as fh:
        _, header, rdr = sniff_reader(fh)
        return header, [dict(zip(header, r)) for r in rdr]


def guess(header, key):
    low = {h.lower(): h for h in header}
    for cand in GUESS[key]:
        if cand in low:
            return low[cand]
    return None


def human(n):
    for u in ["B", "KB", "MB", "GB", "TB"]:
        if n < 1024 or u == "TB":
            return f"{n:.1f} {u}" if u != "B" else f"{n} B"
        n /= 1024


def norm_chrom(c):
    c = str(c).strip()
    return c[3:] if c.lower().startswith("chr") else c


def chrom_sort_key(c):
    return (0, int(c)) if c.isdigit() else (1, {"X": 0, "Y": 1, "M": 2, "MT": 2}.get(c, 9), c)


def safe_name(g):
    return re.sub(r"[^A-Za-z0-9_.-]", "_", g)


def num(x):
    """Turn a CSV string into int/float when possible (keeps JSON small and sortable).

    Genomic coordinates must stay exact. Do not apply significant-figure
    rounding to values >= 1e6; that turns 108440321 into 108400000.
    """
    if x is None:
        return None
    s = str(x).strip()
    if s in ("", "NA", "NaN", "nan", "."):
        return None
    try:
        return int(s)
    except ValueError:
        pass
    try:
        f = float(s)
    except ValueError:
        return s
    if f.is_integer() and abs(f) <= 2**53:
        return int(f)
    if abs(f) >= 1e6:
        return f
    if abs(f) < 1e-3:
        return float(f"{f:.6g}")
    return round(f, 4)


def find_gene_dirs(root, results_dir):
    """Return {gene: folder}. The gene folders are the subfolders of whichever
    directory (outside results/) has the most subfolders that contain files."""
    best, best_n = None, 0
    for dirpath, dirnames, _ in os.walk(root):
        if os.path.abspath(dirpath).startswith(os.path.abspath(results_dir)):
            continue
        if len(dirnames) > best_n:
            best, best_n = dirpath, len(dirnames)
    if not best:
        return None, {}
    genes = {d: os.path.join(best, d) for d in sorted(os.listdir(best))
             if os.path.isdir(os.path.join(best, d))}
    return best, genes


# ---------------------------------------------------------------- 1. inventory
def cmd_inventory(a):
    root = a.root
    results = a.results or os.path.join(root, "results")
    out = []
    p = out.append
    p(f"# Inventory of `{root}`  ({date.today()})\n")

    # results tables
    p("## Summary tables (`results/`)\n")
    p("| file | size | rows | columns |\n|---|---|---|---|")
    tables = {}
    if os.path.isdir(results):
        for f in sorted(os.listdir(results)):
            fp = os.path.join(results, f)
            if not os.path.isfile(fp):
                continue
            size = os.path.getsize(fp)
            if re.search(r"\.(csv|tsv|txt)(\.gz)?$", f):
                with smart_open(fp) as fh:
                    _, header, rdr = sniff_reader(fh)
                    n = sum(1 for _ in rdr)
                tables[f] = header
                p(f"| {f} | {human(size)} | {n:,} | {', '.join(header)} |")
            else:
                p(f"| {f} | {human(size)} | | |")
    else:
        p(f"_No results folder at {results}_")

    # gene folders
    gdir, genes = (a.genes_dir, None) if a.genes_dir else find_gene_dirs(root, results)
    if a.genes_dir:
        genes = {d: os.path.join(gdir, d) for d in sorted(os.listdir(gdir)) if os.path.isdir(os.path.join(gdir, d))}
    p(f"\n## Per-gene folders\n\nFolder: `{gdir}` · {len(genes):,} gene folders\n")
    ext_n, ext_b, patt, patt_example = Counter(), Counter(), Counter(), {}
    per_gene_bytes = Counter()
    for g, d in genes.items():
        for dp, _, fs in os.walk(d):
            for f in fs:
                fp = os.path.join(dp, f)
                b = os.path.getsize(fp)
                ext = f.split(".", 1)[1].lower() if "." in f else "(none)"
                ext_n[ext] += 1; ext_b[ext] += b; per_gene_bytes[g] += b
                rel = os.path.relpath(fp, d).replace(g, "{GENE}")
                rel = re.sub(r"\d+", "#", rel)
                patt[rel] += 1
                patt_example.setdefault(rel, fp)
    p("| extension | files | size |\n|---|---|---|")
    for e, n in ext_n.most_common():
        p(f"| {e} | {n:,} | {human(ext_b[e])} |")
    p(f"\nTotal: {sum(ext_n.values()):,} files, {human(sum(ext_b.values()))}\n")
    p("### File name patterns (gene name → {GENE}, digits → #)\n")
    p("| pattern | files | example header (CSV) |\n|---|---|---|")
    for pt, n in patt.most_common(40):
        hdr = ""
        ex = patt_example[pt]
        if re.search(r"\.(csv|tsv)(\.gz)?$", ex):
            try:
                with smart_open(ex) as fh:
                    hdr = ", ".join(sniff_reader(fh)[1])
            except Exception as e:
                hdr = f"(unreadable: {e})"
        p(f"| `{pt}` | {n:,} | {hdr} |")
    p("\n### Largest genes on disk\n")
    for g, b in per_gene_bytes.most_common(10):
        p(f"- {g}: {human(b)}")

    # reconcile gene sets
    summ_file = next((f for f in tables if "summary" in f.lower()), None)
    hot_file = next((f for f in tables if "hotspot" in f.lower() and f != summ_file), None)
    p("\n## Gene counts\n")
    sets = {"gene folders": set(genes)}
    for label, f in [("summary table", summ_file), ("hotspot table", hot_file)]:
        if f:
            gc = guess(tables[f], "gene")
            if gc:
                _, rows = read_table(os.path.join(results, f))
                sets[f"{label} ({f})"] = {r[gc] for r in rows}
    for k, s in sets.items():
        p(f"- {k}: {len(s):,} genes")
    if len(sets) > 1:
        base = sets["gene folders"]
        for k, s in sets.items():
            if k == "gene folders":
                continue
            only_f, only_t = sorted(base - s), sorted(s - base)
            p(f"- folders not in {k}: {len(only_f):,}" + (f" (e.g. {', '.join(only_f[:15])})" if only_f else ""))
            p(f"- {k} genes without a folder: {len(only_t):,}" + (f" (e.g. {', '.join(only_t[:15])})" if only_t else ""))

    # config skeleton
    var_patt = next((pt for pt, _ in patt.most_common() if re.search(r"variant|score|pred", pt, re.I)
                     and re.search(r"\.(csv|tsv)", pt)), None)
    var_hdr = []
    if var_patt:
        with smart_open(patt_example[var_patt]) as fh:
            var_hdr = sniff_reader(fh)[1]
    cfg = {
        "root": root,
        "results_dir": results,
        "genes_dir": gdir,
        "summary_csv": os.path.join(results, summ_file) if summ_file else "FILL ME",
        "hotspots_csv": os.path.join(results, hot_file) if hot_file else "FILL ME",
        "summary_columns": {"gene": guess(tables.get(summ_file, []), "gene"),
                            "chrom": guess(tables.get(summ_file, []), "chrom"),
                            "show": tables.get(summ_file, [])},
        "hotspot_columns": {k: guess(tables.get(hot_file, []), k) for k in
                            ["gene", "chrom", "start", "end", "strand", "intron", "intron_start", "intron_end",
                             "length", "q", "n_spliceogenic", "max_score"]},
        "variants": {
            "file_glob": (var_patt or "FILL ME e.g. {GENE}_variants.csv").replace("#", "*"),
            "pos": guess(var_hdr, "pos"), "ref": guess(var_hdr, "ref"), "alt": guess(var_hdr, "alt"),
            "score": guess(var_hdr, "score"),
            "min_score": 0.2,
            "max_per_gene": 3000,
            "_note": "web keeps variants with score >= min_score (top max_per_gene by score); 'null' pos disables variants",
        },
        "gtf": None,
        "_gtf_note": "optional: GENCODE/Ensembl GTF (.gz ok) of the same assembly; adds exon coordinates for the intron view",
        "release": {"title": "Genome-wide splice hotspot atlas of disease-associated genes",
                    "version": "v1.0", "assembly": "GRCh38",
                    "related_doi": "10.5281/zenodo.19854065",
                    "web_url": "https://<user>.github.io/<repo>/"},
    }
    with open(a.config_out, "w") as fh:
        json.dump(cfg, fh, indent=2)
    p(f"\n## Next step\n\n`{a.config_out}` was written with guessed column names. Check every value "
      "(a `null` means no guess), then run `web` and `zenodo`.")
    with open(a.out, "w") as fh:
        fh.write("\n".join(out) + "\n")
    print("\n".join(out))
    print(f"\nwrote {a.out} and {a.config_out}")


# ---------------------------------------------------------------- shared loading
def load_cfg(path):
    with open(path) as fh:
        cfg = json.load(fh)
    for k in ["summary_csv", "hotspots_csv", "genes_dir"]:
        if not cfg.get(k) or "FILL" in str(cfg[k]):
            sys.exit(f"config: set {k}")
    return cfg


def load_tables(cfg):
    sc, hc = cfg["summary_columns"], cfg["hotspot_columns"]
    _, srows = read_table(cfg["summary_csv"])
    _, hrows = read_table(cfg["hotspots_csv"])
    for need in ["gene"]:
        if not sc.get(need) or not hc.get(need):
            sys.exit(f"config: summary_columns.{need} and hotspot_columns.{need} are required")
    summary = {r[sc["gene"]]: r for r in srows}
    hot = defaultdict(list)
    for r in hrows:
        hot[r[hc["gene"]]].append(r)
    return summary, hot


def load_intron_interior(cfg):
    """Map gene -> intron_name -> n_spliceogenic_interior for concentration."""
    path = cfg.get("intron_summary_csv")
    if not path or not os.path.isfile(path):
        return {}
    gene_col = cfg["hotspot_columns"].get("gene") or "gene"
    intron_col = cfg["hotspot_columns"].get("intron") or "intron_name"
    out = defaultdict(dict)
    _, rows = read_table(path)
    for r in rows:
        g = r.get(gene_col)
        intron = r.get(intron_col)
        if not g or not intron:
            continue
        n = num(r.get("n_spliceogenic_interior"))
        if n is None:
            continue
        out[g][str(intron)] = int(n)
    return out


def gene_chrom(g, cfg, summary, hot):
    sc, hc = cfg["summary_columns"], cfg["hotspot_columns"]
    if sc.get("chrom") and g in summary and summary[g].get(sc["chrom"]):
        return norm_chrom(summary[g][sc["chrom"]])
    if hc.get("chrom") and hot.get(g):
        return norm_chrom(hot[g][0][hc["chrom"]])
    return "unplaced"


ATTR = re.compile(r'(\S+) "([^"]*)"')


def read_gtf(path, wanted):
    """Optional fallback. Prefer load_scored_models() for this screen."""
    tx = {}
    with smart_open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "exon":
                continue
            at = defaultdict(list)
            for k, v in ATTR.findall(f[8]):
                at[k].append(v)
            g = (at.get("gene_name") or at.get("gene_id"))[0]
            if g not in wanted:
                continue
            t = tx.setdefault(at["transcript_id"][0], {"gene": g, "chrom": norm_chrom(f[0]), "strand": f[6],
                                                        "exons": [], "tags": set(at.get("tag", []))})
            t["exons"].append([int(f[3]), int(f[4])])
    best = {}
    rank = lambda t: ("MANE_Select" in t["tags"], "Ensembl_canonical" in t["tags"],
                      sum(e - s + 1 for s, e in t["exons"]))
    for tid, t in tx.items():
        t["id"] = tid
        if t["gene"] not in best or rank(t) > rank(best[t["gene"]]):
            best[t["gene"]] = t
    for t in best.values():
        t["exons"].sort()
        t.pop("tags")
    return best


def parse_intron_number(name):
    if name is None:
        return None
    m = INTRON_NAME_RE.match(str(name).strip())
    return int(m.group(1)) if m else num(name)


def load_transcript_meta(gdir):
    path = os.path.join(gdir, "transcript.json")
    if not os.path.isfile(path):
        return {}
    with open(path) as fh:
        return json.load(fh)


def load_intron_bounds(gdir):
    """Map intron_name to 1-based inclusive bounds from *_introns.tsv (0-based half-open)."""
    matches = sorted(glob.glob(os.path.join(gdir, "*_introns.tsv")))
    if not matches:
        return {}
    out = {}
    with open(matches[0], newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            name = row.get("intron_name") or row.get("name")
            if not name:
                continue
            try:
                start0 = int(row["start_hg38"])
                end0 = int(row["end_hg38"])
            except (KeyError, ValueError, TypeError):
                continue
            out[name] = {
                "intron_start": start0 + 1,
                "intron_end": end0,
                "strand": row.get("strand"),
                "length": int(row["length_bp"]) if row.get("length_bp") not in (None, "") else end0 - start0,
            }
    return out


def load_scored_models(genes_dir, genes, refseq_path=None):
    """Build transcript models from transcript.json and optional UCSC ncbiRefSeq exons."""
    models = {}
    wanted_tx = {}
    for g in genes:
        gdir = os.path.join(genes_dir, g)
        meta = load_transcript_meta(gdir) if os.path.isdir(gdir) else {}
        if not meta:
            continue
        tid = meta.get("transcript") or meta.get("transcript_id")
        models[g] = {
            "id": tid,
            "chrom": norm_chrom(meta.get("chrom", "")),
            "strand": meta.get("strand"),
            "exons": [],
            "exon_count": meta.get("exon_count"),
            "status": meta.get("status"),
        }
        if tid:
            wanted_tx[tid] = g
    if not refseq_path or not os.path.isfile(refseq_path) or not wanted_tx:
        return models
    # Prefer primary-assembly rows. Duplicate transcript ids (alt/fix patches)
    # must not inflate the "found" count or the scan stops before chrX/Y genes.
    pending = dict(wanted_tx)
    with open(refseq_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 11:
                continue
            tid = f[1]
            if tid not in pending:
                continue
            chrom = f[2]
            g = pending[tid]
            starts = [int(x) for x in f[9].rstrip(",").split(",") if x]
            ends = [int(x) for x in f[10].rstrip(",").split(",") if x]
            # UCSC exon coords are 0-based half-open; site uses 1-based inclusive.
            exons = [[s + 1, e] for s, e in zip(starts, ends)]
            models[g]["exons"] = exons
            models[g]["chrom"] = norm_chrom(chrom)
            models[g]["strand"] = f[3]
            # Keep the first primary-assembly hit; allow a later primary to
            # replace an earlier alt/fix row for the same transcript.
            if "_" not in chrom:
                pending.pop(tid)
            if not pending:
                break
    return models


# ---------------------------------------------------------------- 2. web
def gene_variant_files(gdir, gene, pattern):
    pat = pattern.replace("{GENE}", glob.escape(gene))
    hits = []
    for dp, _, fs in os.walk(gdir):
        for f in fs:
            rel = os.path.relpath(os.path.join(dp, f), gdir).replace(os.sep, "/")
            if fnmatch.fnmatchcase(rel, pat):
                hits.append(os.path.join(dp, f))
    return hits


def enrich_intervals(hrows, hc, bounds, interior_by_intron):
    """Build interval records from shortlist rows. Concentration uses intron interior."""
    intervals = []
    for r in hrows:
        iv = {k: num(r.get(v)) for k, v in hc.items() if v}
        raw_intron = r.get(hc["intron"]) if hc.get("intron") else None
        intron_n = parse_intron_number(raw_intron)
        if intron_n is not None:
            iv["intron"] = intron_n
        b = bounds.get(str(raw_intron)) if raw_intron is not None else None
        if b:
            iv["intron_start"] = b["intron_start"]
            iv["intron_end"] = b["intron_end"]
            if b.get("strand"):
                iv["strand"] = b["strand"]
        n_spl = iv.get("n_spliceogenic")
        n_int = interior_by_intron.get(str(raw_intron)) if raw_intron is not None else None
        if n_int is None and intron_n is not None:
            n_int = interior_by_intron.get(f"intron_{intron_n}")
        if n_spl is not None and n_int:
            iv["concentration"] = round(float(n_spl) / float(n_int), 4)
            iv["n_spliceogenic_interior"] = int(n_int)
        elif n_int is not None:
            iv["n_spliceogenic_interior"] = int(n_int)
            iv["concentration"] = None
        else:
            iv["concentration"] = None
        intervals.append(iv)
    return intervals


def finalize_intervals(intervals):
    """Sort by min_BH_q, mark default and top-in-intron (best q per intron)."""
    def sort_key(iv):
        q = iv.get("q")
        n = iv.get("n_spliceogenic")
        n = -int(n) if n is not None else 0
        return (q is None, q if q is not None else 0, n)

    intervals = sorted(intervals, key=sort_key)
    seen = set()
    for i, iv in enumerate(intervals):
        intron = iv.get("intron")
        iv["top_in_intron"] = intron not in seen and intron is not None
        if intron is not None:
            seen.add(intron)
        iv["default"] = i == 0
        iv["id"] = f"IV-{i + 1}"
    return intervals


def build_gene(job):
    g, gdir, srow, hrows, model, interior_by_intron, cfg = job
    hc, vc = cfg["hotspot_columns"], cfg["variants"]
    bounds = load_intron_bounds(gdir) if gdir else {}
    intervals = enrich_intervals(hrows, hc, bounds, interior_by_intron or {})
    variants, scanned = [], 0
    if vc.get("pos") and vc.get("score") and gdir:
        for fp in gene_variant_files(gdir, g, vc["file_glob"]):
            with smart_open(fp) as fh:
                _, header, rdr = sniff_reader(fh)
                ix = {h: i for i, h in enumerate(header)}
                if vc["pos"] not in ix or vc["score"] not in ix:
                    continue
                pi, si = ix[vc["pos"]], ix[vc["score"]]
                ri = ix.get(vc.get("ref") or "")
                ai = ix.get(vc.get("alt") or "")
                for row in rdr:
                    scanned += 1
                    try:
                        s = float(row[si])
                        pos = int(float(row[pi]))
                    except (ValueError, IndexError):
                        continue
                    if s < vc["min_score"]:
                        continue
                    v = [pos, round(s, 3)]
                    if ri is not None and ai is not None:
                        v += [row[ri], row[ai]]
                    variants.append(v)
                    for iv in intervals:
                        start, end = iv.get("start"), iv.get("end")
                        if start is None or end is None:
                            continue
                        if start <= pos <= end:
                            prev = iv.get("max_score")
                            iv["max_score"] = s if prev is None else max(prev, s)
        variants.sort(key=lambda v: -v[1])
        variants = sorted(variants[: vc["max_per_gene"]], key=lambda v: v[0])
        for iv in intervals:
            if iv.get("max_score") is not None:
                iv["max_score"] = round(float(iv["max_score"]), 3)
    intervals = finalize_intervals(intervals)
    n_top = sum(1 for iv in intervals if iv.get("top_in_intron"))
    doc = {"gene": g,
           "summary": {k: num(v) for k, v in (srow or {}).items() if k in (cfg["summary_columns"].get("show") or [])},
           "intervals": intervals, "variants": variants,
           "variant_columns": ["pos", "score", "ref", "alt"][: len(variants[0])] if variants else [],
           "interval_columns": [
               "id", "chrom", "start", "end", "intron", "length", "q",
               "n_spliceogenic", "concentration", "max_score", "top_in_intron",
               "default", "intron_start", "intron_end", "strand",
           ]}
    if model:
        if model.get("id"):
            doc["transcript"] = model["id"]
        if model.get("chrom"):
            doc["chrom"] = model["chrom"]
        if model.get("strand"):
            doc["strand"] = model["strand"]
        if model.get("exons"):
            doc["exons"] = model["exons"]
    elif gdir:
        meta = load_transcript_meta(gdir)
        if meta.get("transcript"):
            doc["transcript"] = meta["transcript"]
        if meta.get("chrom"):
            doc["chrom"] = norm_chrom(meta["chrom"])
        if meta.get("strand"):
            doc["strand"] = meta["strand"]
    out = cfg["_out"]
    sname = safe_name(g)
    with open(os.path.join(out, "genes", f"{sname}.json"), "w") as fh:
        json.dump(doc, fh, separators=(",", ":"))
    dl_cols = [
        "id", "chrom", "start", "end", "intron", "length", "q",
        "n_spliceogenic", "n_spliceogenic_interior", "concentration",
        "max_score", "top_in_intron", "intron_start", "intron_end", "strand",
    ]
    with open(os.path.join(out, "downloads", f"{sname}_intervals.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(dl_cols)
        for iv in intervals:
            w.writerow([iv.get(c, "") for c in dl_cols])
    # Remove stale hotspot download name from earlier builds.
    stale = os.path.join(out, "downloads", f"{sname}_hotspots.csv")
    if os.path.isfile(stale):
        os.unlink(stale)
    return g, sname, len(intervals), n_top, len(variants), scanned


def cmd_web(a):
    cfg = load_cfg(a.config)
    out = a.out
    os.makedirs(os.path.join(out, "genes"), exist_ok=True)
    os.makedirs(os.path.join(out, "downloads"), exist_ok=True)
    cfg["_out"] = out
    summary, hot = load_tables(cfg)
    interior = load_intron_interior(cfg)
    gene_dirs = {d: os.path.join(cfg["genes_dir"], d) for d in os.listdir(cfg["genes_dir"])
                 if os.path.isdir(os.path.join(cfg["genes_dir"], d))}
    # Publish only genes that were scored into the summary table.
    genes = sorted(summary)
    if cfg.get("gtf"):
        models = read_gtf(cfg["gtf"], set(genes))
        print(f"web: loaded {len(models):,} transcript models from GTF {cfg['gtf']}")
    else:
        models = load_scored_models(cfg["genes_dir"], genes, cfg.get("refseq"))
        n_exons = sum(1 for m in models.values() if m.get("exons"))
        print(f"web: {len(models):,} scored transcripts from gene folders; "
              f"{n_exons:,} with exons from refseq")
    print(f"web: {len(genes):,} genes in summary, {sum(1 for g in genes if hot.get(g)):,} with intervals; "
          f"interior stats for {len(interior):,} genes")
    jobs = [
        (g, gene_dirs.get(g), summary.get(g), hot.get(g, []), models.get(g), interior.get(g, {}), cfg)
        for g in genes
    ]
    index, done = {}, 0
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        futs = [ex.submit(build_gene, j) for j in jobs]
        for f in as_completed(futs):
            g, sname, n_iv, n_top, nv, scanned = f.result()
            done += 1
            index[g] = {
                "f": sname,
                "chrom": gene_chrom(g, cfg, summary, hot),
                "h": n_iv,
                "t": n_top,
                "v": nv,
            }
            if done % 250 == 0 or done == len(jobs):
                print(f"  {done:,}/{len(jobs):,} genes")
    sc = cfg["summary_columns"]
    for g in index:   # add a few summary fields to the index for search results
        row = summary.get(g, {})
        index[g]["s"] = {k: num(row.get(k)) for k in (sc.get("show") or [])[:6] if k != sc.get("gene")}
    rel = cfg.get("release", {})
    with open(os.path.join(out, "index.json"), "w") as fh:
        json.dump({"meta": {"built": str(date.today()), "assembly": rel.get("assembly"),
                            "version": rel.get("version"), "min_variant_score": cfg["variants"]["min_score"],
                            "n_genes": len(index),
                            "interval_unit": "interior significant interval",
                            "interval_sort": "min_BH_q",
                            "top_interval": "best min_BH_q per intron"},
                   "genes": index}, fh, separators=(",", ":"))
    # size report
    sizes = []
    for dp, _, fs in os.walk(out):
        for f in fs:
            fp = os.path.join(dp, f)
            sizes.append((os.path.getsize(fp), os.path.relpath(fp, out)))
    total = sum(s for s, _ in sizes)
    print(f"\nwrote {len(sizes):,} files, {human(total)} to {out}/")
    print("largest:", ", ".join(f"{p} ({human(s)})" for s, p in sorted(sizes, reverse=True)[:5]))
    if total > 900 * 1024**2:
        print("WARNING: over ~900 MB. GitHub Pages sites are limited to 1 GB; raise variants.min_score "
              "or lower variants.max_per_gene.")
    if sizes and max(s for s, _ in sizes) > 50 * 1024**2:
        print("WARNING: a file is over 50 MB (GitHub warns at 50 MB, rejects at 100 MB).")


# ---------------------------------------------------------------- 3. zenodo
STORE_EXT = {".png", ".jpg", ".jpeg", ".gz", ".zip", ".bam", ".bai", ".pdf"}


def add_file(zf, fp, arc, manifest, zipname):
    ext = os.path.splitext(fp)[1].lower()
    comp = zipfile.ZIP_STORED if ext in STORE_EXT else zipfile.ZIP_DEFLATED
    info = zipfile.ZipInfo.from_file(fp, arc)
    info.compress_type = comp
    h, n = hashlib.sha256(), 0
    with open(fp, "rb") as src, zf.open(info, "w", force_zip64=True) as dst:
        for chunk in iter(lambda: src.read(1 << 20), b""):
            h.update(chunk); dst.write(chunk); n += len(chunk)
    manifest.append((zipname, arc, n, h.hexdigest()))


def md5(fp):
    h = hashlib.md5()
    with open(fp, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def classify_unscored(gdir):
    """Return a short reason for a gene folder that is absent from the summary."""
    meta = load_transcript_meta(gdir)
    status = meta.get("status")
    if status == "single_exon":
        return "single_exon"
    if status == "ok" and glob.glob(os.path.join(gdir, "*_introns.tsv")):
        return "no_atlas_rows"
    if status:
        return str(status)
    if not meta:
        return "no_transcript_json"
    return "unscored"


def write_not_scored(path, genes_dir, summary_genes):
    rows = []
    for name in sorted(os.listdir(genes_dir)):
        gdir = os.path.join(genes_dir, name)
        if not os.path.isdir(gdir) or name in summary_genes:
            continue
        meta = load_transcript_meta(gdir)
        rows.append({
            "gene": name,
            "reason": classify_unscored(gdir),
            "transcript": meta.get("transcript", ""),
            "chrom": meta.get("chrom", ""),
            "status": meta.get("status", ""),
            "exon_count": meta.get("exon_count", ""),
            "intron_count": meta.get("intron_count", ""),
        })
    with open(path, "w", newline="") as fh:
        fields = ["gene", "reason", "transcript", "chrom", "status", "exon_count", "intron_count"]
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    return rows


def cmd_zenodo(a):
    cfg = load_cfg(a.config)
    out = a.out
    os.makedirs(out, exist_ok=True)
    summary, hot = load_tables(cfg)
    rel = cfg.get("release", {})
    ver = rel.get("version", "v1")
    gdir = cfg["genes_dir"]
    # Only scored genes enter chromosome archives.
    genes = sorted(g for g in summary if os.path.isdir(os.path.join(gdir, g)))
    missing_dirs = sorted(set(summary) - set(genes))
    if missing_dirs:
        print(f"WARNING: {len(missing_dirs)} summary genes have no folder "
              f"(e.g. {', '.join(missing_dirs[:10])})")
    skip_ext = set() if a.include_svg else {".svg"}
    manifest, zips = [], []
    results_files = cfg.get("results_files") or list(DEFAULT_RESULTS_FILES)

    # summary tables: allowlisted CSVs only (skip logs)
    zname = f"summary_tables_{ver}.zip"
    with zipfile.ZipFile(os.path.join(out, zname), "w", zipfile.ZIP_DEFLATED, compresslevel=6) as zf:
        for f in results_files:
            fp = os.path.join(cfg["results_dir"], f)
            if not os.path.isfile(fp):
                print(f"WARNING: missing results file {fp}")
                continue
            add_file(zf, fp, f"results/{f}", manifest, zname)
    zips.append(zname)
    print(f"  {zname}")

    # per-chromosome gene archives
    by_chr = defaultdict(list)
    for g in genes:
        by_chr[gene_chrom(g, cfg, summary, hot)].append(g)
    for c in sorted(by_chr, key=chrom_sort_key):
        zname = f"genes_chr{c}_{ver}.zip" if c != "unplaced" else f"genes_unplaced_{ver}.zip"
        with zipfile.ZipFile(os.path.join(out, zname), "w", zipfile.ZIP_DEFLATED, compresslevel=6) as zf:
            for g in by_chr[c]:
                for dp, _, fs in os.walk(os.path.join(gdir, g)):
                    for f in sorted(fs):
                        if os.path.splitext(f)[1].lower() in skip_ext:
                            continue
                        fp = os.path.join(dp, f)
                        arc = f"genes/chr{c}/{g}/" + os.path.relpath(fp, os.path.join(gdir, g)).replace(os.sep, "/")
                        add_file(zf, fp, arc, manifest, zname)
        zips.append(zname)
        print(f"  {zname}: {len(by_chr[c]):,} genes, {human(os.path.getsize(os.path.join(out, zname)))}")

    not_scored = write_not_scored(os.path.join(out, "not_scored.tsv"), gdir, set(summary))
    print(f"  not_scored.tsv: {len(not_scored):,} genes")
    reason_counts = Counter(r["reason"] for r in not_scored)
    for reason, n in reason_counts.most_common():
        print(f"    {reason}: {n:,}")

    with open(os.path.join(out, "MANIFEST.tsv"), "w") as fh:
        fh.write("archive\tpath\tbytes\tsha256\n")
        for r in manifest:
            fh.write("\t".join(map(str, r)) + "\n")
    with open(os.path.join(out, "gene_to_archive.tsv"), "w") as fh:
        fh.write("gene\tchrom\tarchive\n")
        for c, gs in by_chr.items():
            z = f"genes_chr{c}_{ver}.zip" if c != "unplaced" else f"genes_unplaced_{ver}.zip"
            for g in sorted(gs):
                fh.write(f"{g}\t{c}\t{z}\n")
    with open(os.path.join(out, "checksums.md5"), "w") as fh:
        for z in zips:
            fh.write(f"{md5(os.path.join(out, z))}  {z}\n")
    write_readme(out, cfg, zips, len(genes), manifest, a.include_svg, reason_counts)

    files = sorted(os.listdir(out))
    total = sum(os.path.getsize(os.path.join(out, f)) for f in files)
    print(f"\nrelease: {len(files)} files, {human(total)} in {out}/")
    if len(files) > 100:
        print("WARNING: Zenodo allows at most 100 files per record.")
    if total > 50e9:
        print("WARNING: over Zenodo's 50 GB default. Drop SVGs (omit --include-svg) or request a quota increase.")


def write_readme(out, cfg, zips, n_genes, manifest, svg, reason_counts=None):
    rel = cfg.get("release", {})
    by_ext = Counter(os.path.splitext(p)[1].lower() for _, p, _, _ in manifest)
    reason_counts = reason_counts or Counter()
    not_scored_line = ", ".join(f"{k} {v:,}" for k, v in reason_counts.most_common()) or "none"
    lines = [
        f"# {rel.get('title', 'Splice hotspot atlas')} ({rel.get('version', '')})", "",
        f"Genome assembly: {rel.get('assembly', '')}. Scored genes: {n_genes:,}. Built {date.today()}.", "",
        f"Interactive browser: {rel.get('web_url', '')}",
        f"Related record (RNA-seq BAMs): https://doi.org/{rel.get('related_doi', '')}", "",
        "## Files", "",
        "| file | contents |", "|---|---|",
        f"| summary_tables_{rel.get('version','')}.zip | genome-wide summary and hotspot tables (`results/`) |",
        "| genes_chr*_*.zip | per-gene tables and plots for scored genes: `genes/chr<N>/<GENE>/` |",
        "| gene_to_archive.tsv | which archive holds each scored gene |",
        "| not_scored.tsv | disease-gene folders omitted from scoring (reason per gene) |",
        "| MANIFEST.tsv | every file in every archive with size and SHA-256 |",
        "| checksums.md5 | MD5 of each archive (compare with Zenodo's listed checksums) |", "",
        "File types in the archives: " + ", ".join(f"{e or '(none)'} {n:,}" for e, n in by_ext.most_common()) + ".",
        "" if svg else "Vector (SVG) plots are not included in this release; PNGs are.", "",
        "## Genes not scored", "",
        f"Counts by reason: {not_scored_line}.",
        "single_exon: chosen RefSeq transcript has no intron.",
        "no_atlas_rows: introns were extracted but the Atlas returned no SNVs for those intervals.", "",
        "## Column descriptions", "", "[FILL IN: one line per column of each summary table]", "",
        "## Methods", "", "[FILL IN: prediction model and version, score threshold, hotspot test, BH correction]", "",
        "## How to cite", "", "[FILL IN after the DOI is reserved]", "",
        "## License", "", "[e.g. CC BY 4.0]", "",
    ]
    with open(os.path.join(out, "README.md"), "w") as fh:
        fh.write("\n".join(lines))


# ---------------------------------------------------------------- main
def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    i = sub.add_parser("inventory", help="inspect folders, write inventory.md + config.json")
    i.add_argument("--root", required=True, help="top folder, e.g. genome_wide")
    i.add_argument("--results", help="summary table folder (default <root>/results)")
    i.add_argument("--genes-dir", help="folder holding one subfolder per gene (default: auto-detect)")
    i.add_argument("--out", default="inventory.md")
    i.add_argument("--config-out", default="config.json")
    w = sub.add_parser("web", help="build website data")
    w.add_argument("--config", default="config.json")
    w.add_argument("--out", default="docs/data")
    w.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    z = sub.add_parser("zenodo", help="build Zenodo release archives")
    z.add_argument("--config", default="config.json")
    z.add_argument("--out", default="release")
    z.add_argument("--include-svg", action="store_true", help="also pack SVG plots (bigger release)")
    a = p.parse_args()
    {"inventory": cmd_inventory, "web": cmd_web, "zenodo": cmd_zenodo}[a.cmd](a)


if __name__ == "__main__":
    main()
