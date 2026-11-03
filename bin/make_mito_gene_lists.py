#!/usr/bin/env python
"""Build assets/mito_genes_{hs,mm}.tsv: every gene on the mitochondrial chromosome, from Ensembl and NCBI Gene.

Run: python bin/make_mito_gene_lists.py [--outdir assets]
The mitochondrial genome encodes 37 genes in both species: 13 protein-coding, 2 rRNA and 22 tRNA genes.
Ensembl REST (/overlap/region/<species>/MT) gives the Ensembl IDs, nomenclature symbols (MT-CO1, mt-Co1)
and biotypes. NCBI Gene gene_info gives the NCBI ID and NCBI Symbol (COX1), which Gemma imports; rows are
joined on the nomenclature symbol (Symbol_from_nomenclature_authority), which also drops NCBI's extra
chrM records for other mitogenomes filed under the same taxon. All 37 genes are listed. Most 10x/Gemma
gene sets carry only the 13 protein-coding and 2 rRNA genes, and the pipeline matches whatever is present.
Columns: ensembl_id, symbol, ncbi_gene_id, gemma_symbol, gene_type. A sidecar mito_genes_<org>.source.json
records the Ensembl release, URLs, NCBI Last-Modified header, sha256 of each download and the download time.
Four requests in all, one second apart.
"""
import argparse, csv, datetime, gzip, hashlib, io, json, os, sys, time, urllib.request

ORGS = {"hs": ("human", "Homo_sapiens", 9606, 16569), "mm": ("mouse", "Mus_musculus", 10090, 16299)}
BIOTYPE = {"protein_coding": "protein-coding", "Mt_rRNA": "rRNA", "Mt_tRNA": "tRNA"}
EXPECTED = {"protein-coding": 13, "rRNA": 2, "tRNA": 22}
UA = {"User-Agent": "sc-annotation-pipeline make_mito_gene_lists.py", "Content-Type": "application/json"}

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--outdir", default=os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "assets"))
a = ap.parse_args()


def get(url, tries=3):
    for i in range(tries):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=300) as r:
                body = r.read()
            if body:
                time.sleep(1)
                return body, r.headers
        except OSError as e:
            err = e
        time.sleep(5 * (i + 1))
    sys.exit(f"failed to fetch {url}")


release = json.loads(get("https://rest.ensembl.org/info/data")[0])["releases"][0]
for org, (species, name, taxid, mt_len) in ORGS.items():
    ens_url = f"https://rest.ensembl.org/overlap/region/{species}/MT:1-{mt_len}?feature=gene"
    ens_raw, _ = get(ens_url)
    ens = [g for g in json.loads(ens_raw) if g["seq_region_name"] == "MT"]
    ncbi_url = f"https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/{name}.gene_info.gz"
    ncbi_raw, hdr = get(ncbi_url)
    ncbi = {}
    for d in csv.DictReader(io.TextIOWrapper(gzip.GzipFile(fileobj=io.BytesIO(ncbi_raw)), "utf-8"), delimiter="\t"):
        if d["#tax_id"] == str(taxid) and d["chromosome"] == "MT" and d["Symbol_from_nomenclature_authority"] != "-":
            ncbi[d["Symbol_from_nomenclature_authority"]] = d
    rows = []
    for g in ens:
        sym, n = g["external_name"], ncbi.get(g["external_name"])
        if n is None:
            sys.exit(f"{org}: Ensembl gene {g['id']} {sym} has no NCBI Gene record with that nomenclature symbol")
        rows.append({"ensembl_id": g["id"], "symbol": sym, "ncbi_gene_id": n["GeneID"], "gemma_symbol": n["Symbol"],
                     "gene_type": BIOTYPE.get(g["biotype"], g["biotype"])})
    counts = {t: sum(r["gene_type"] == t for r in rows) for t in EXPECTED}
    if counts != EXPECTED or len(rows) != sum(EXPECTED.values()):
        sys.exit(f"{org}: unexpected chrM gene set {counts}, {len(rows)} genes (expected {EXPECTED})")
    rows.sort(key=lambda r: (list(EXPECTED).index(r["gene_type"]), r["symbol"]))
    out = os.path.join(a.outdir, f"mito_genes_{org}.tsv")
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(rows)
    with open(os.path.join(a.outdir, f"mito_genes_{org}.source.json"), "w") as f:
        json.dump({"ensembl_release": release, "ensembl_url": ens_url, "ensembl_sha256": hashlib.sha256(ens_raw).hexdigest(),
                   "ncbi_url": ncbi_url, "ncbi_last_modified": hdr.get("Last-Modified"),
                   "ncbi_sha256": hashlib.sha256(ncbi_raw).hexdigest(),
                   "downloaded": datetime.datetime.now().isoformat(timespec="seconds"), "by_type": counts}, f, indent=1)
    print(f"{org}: {len(rows)} chrM genes {counts} -> {out}")
