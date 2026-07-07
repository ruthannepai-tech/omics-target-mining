# omics-target-mining

An end-to-end, reproducible pipeline for mining **public omics datasets**
(GEO, ArrayExpress, PRIDE, CELLxGENE) to nominate and prioritize
**therapeutic targets** for a disease — going from *"here is a condition"* to
*"here is a ranked, druggability-annotated target shortlist"* using only
public data.

It is packaged as a [Claude Science](https://claude.com) **skill**
(`SKILL.md` + a `kernel.py` sidecar), but the two core statistics in
`kernel.py` are plain, dependency-light Python functions that work on their
own in any analysis.

## What it does

Ten stages, each ending by saving artifacts:

1. **Dataset inventory** — query GEO / ArrayExpress / PRIDE / CELLxGENE, build one master table, audit for keyword-only false positives.
2. **Tier & prioritize** — score datasets on relevance, analytical role, tractability (Tier 1 bulk case-control → Tier 3 complementary).
3. **Retrieve & QC Tier-1 bulk** — harmonize IDs to gene symbols, normalize, PCA + canonical-marker sanity checks.
4. **Per-study differential expression** — variance-moderated Welch *t*-test with BH-FDR (`moderated_de`).
5. **Cross-study meta-analysis** — DerSimonian–Laird random-effects pooling with heterogeneity stats (`meta_analysis`).
6. **Single-cell validation** — scanpy QC → HVG → Harmony → Leiden → composition shift + per-cell-type DE.
7. **Pathway / regulator / cell-cell communication** — Enrichr (GO_BP, Reactome, Hallmark, TF libraries) + LIANA.
8. **Multi-omics integration** — per-gene evidence table across independently-varying layers.
9. **Druggability** — Open Targets GraphQL: tractability, clinical candidates, subcellular location, target class.
10. **Consolidate** — Tier-A shortlist with mechanism rationale and a suggested modality per target.

## The `kernel.py` helpers (usable standalone)

Two functions, depending only on `numpy`, `pandas`, `scipy`, and
`statsmodels`:

```python
from kernel import moderated_de, meta_analysis

# Stage 4: per-study DE on log2 expression (genes x samples)
de = moderated_de(expr, labels, case="disease", control="control")
#   -> DataFrame: gene, log2FC, t, pval, padj  (variance-moderated Welch t + BH-FDR)

# Stage 5: pool per-study effects (genes x studies of log2FC and SE)
meta = meta_analysis(fc_df, se_df, min_studies=3)
#   -> DataFrame: pooled_log2FC, se, z, pval, padj, Q, I2, k, n_up, n_dn, direction
#      (DerSimonian-Laird random-effects meta-analysis)
```

## Install

```bash
pip install -r requirements.txt
```

The full pipeline additionally uses `requests`, `openpyxl`, `GEOparse`,
`gseapy`, `mygene`, and — for single-cell stages — `scanpy`, `leidenalg`,
`harmonypy`, and `liana`.

## Using it as a Claude Science skill

Drop the `omics-target-mining/` folder (containing `SKILL.md` and
`kernel.py`) into your skills directory, or publish it via the skill SDK.
When the skill loads, `moderated_de` and `meta_analysis` are auto-defined in
the kernel and `SKILL.md` drives the 10-stage workflow.

## License

MIT — see [LICENSE](LICENSE).
