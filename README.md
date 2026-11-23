# Nonsense-mediated decay escape in neurological genes and reported phenotype of truncating variants

**Taner Karagöl**<sup>1,\*,‡</sup> **& Alper Karagöl**<sup>2,\*</sup>

<sup>1</sup> McGovern Institute for Brain Research, Massachusetts Institute of Technology, Cambridge, Massachusetts, United States

<sup>2</sup> Department of Mechanical Engineering, Massachusetts Institute of Technology, Cambridge, Massachusetts, United States

<sup>\*</sup> These authors contributed equally.

<sup>‡</sup> To whom correspondence should be addressed.

## Abstract

Truncating variants are read under a loss-of-function model in which position matters only through escape from nonsense-mediated decay. In single genes, stop position can select which disease follows, and the reading frame can decide pathogenicity. We tested both across the 1,376 genes that ClinVar links to a neurological condition, with panel rules fixed before any result. We took phenotypes from each pathogenic truncating variant's submissions, and gnomAD v4 predicted loss-of-function variants of the same gene formed the population control. Every test was repeated without population variants on another transcript, without variants in low-complexity sequence, and without both. In the 110 genes with both a neurological and another reported phenotype, position has no detectable effect on which is reported. No gene separates in any sub-unit or under any escape-window definition or restriction, and all 440 whole-panel gene-window comparisons return q = 1.000. The reading frame separates pathogenic from population truncations in nine of 667 genes, six of them through population frameshifts in low-complexity sequence; ARID1A, ATM and CIC alone separate under every restriction. The frameshift-to-stop distance separates in five of 436 genes, four with shorter pathogenic tails, in SYNGAP1 alone under both restrictions and in none under gnomAD's low-complexity flag. Escape-window membership separates pathogenic from population truncations in eight genes, six under every restriction. In four of these six, the population has window truncations where few or no pathogenic ones fall. Position therefore bears on whether a late truncation is tolerated and has no detectable effect on which phenotype is reported.

## Keywords

Truncating variants, Nonsense-mediated decay, Frameshift, PVS1, Low-complexity regions, Neurodevelopmental disorders, Neurodegeneration

## Repository layout

```
neurological-truncation-phenotype/
├── README.md
├── LICENSE
├── Supplementary_Information.pdf
├── codes/
│   ├── run_all.sh                     the whole analysis in order, from the dated ClinVar files
│   ├── gate.py                        the forty panel cuts on which the panel rules were fixed
│   ├── grouping_report.py             condition names the gene-condition map attaches to many genes
│   ├── panel_v2.py                    the whole-panel term lists and the cut of the gene-condition map
│   ├── run_panel.py                   cuts a unit's panel with the cached UniProt answers and the grouping rule
│   ├── uniprot_rows.py                queries UniProt once per symbol and writes uniprot_rows.json
│   ├── uniprot_rows.json              the UniProt answers the panel cuts are resolved with
│   ├── area.py                        the panel's test of a reported condition name
│   ├── panel_symbols.py               the gene symbols of the forty panel cuts
│   ├── vs_filter.py                   the variant summary cut to those symbols
│   ├── ss_filter.py                   the submission summary cut to the variants kept
│   ├── gnomad_union.py                the gene list of the gnomAD pull
│   ├── gnomad_pull.py                 gnomAD v4 predicted loss-of-function variants and exon structure per gene
│   ├── gnomad_pull_durable.py         runs gnomad_pull.py with every API answer stored as it arrives
│   ├── run_unit.sh                    the pipeline of one unit
│   ├── review_filter.py               per-gene truncating counts under a submission-quality filter
│   ├── fsd2.py                        frameshift-to-stop distance on the full RefSeq transcript
│   ├── export_variant_table.py        one row per pathogenic and population truncation of the isolated tier
│   ├── frame_dependence.py            the four frame tests per gene
│   ├── fsd_summary.py                 pooled novel-tail statistics
│   ├── fsd_q.py                       corrected q values of the per-gene tail tests
│   ├── dual_phenotypes.py             the genes with a second condition, under the broad and the strict rule
│   ├── pilot_positions.py             position against the reported phenotype, by median position
│   ├── pilot_windows.py               escape-window membership under the four window definitions
│   ├── labels_submission_level.py     the reported phenotype per variant from the submission records
│   ├── position_q.py                  corrected q values of the position axis
│   ├── pooled_exclusions.py           each pooled test with ATM and with the carrying gene removed
│   ├── run_restrictions.sh            the transcript, low-complexity and combined restriction runs
│   ├── gnomad_restrict.py             the gnomAD records under each restriction
│   ├── lcr_intervals.py               membership of a locus in the low-complexity intervals
│   ├── vs_lcr_filter.py               the variant summary outside low-complexity sequence
│   ├── lcr_counts.py                  low-complexity counts among the isolated-tier variants
│   ├── label_audit.py                 which reported conditions fall on each side of the label
│   ├── genomewide_power.py            the pathogenic frameshift share of every ClinVar gene
│   ├── genomewide_rank.py             where the frame genes sit in that distribution
│   ├── pathogenic_accessions.py       the RefSeq transcript on which ClinVar names each pathogenic truncation
│   ├── collate.py                     one summary per unit and restriction, results_summary.json
│   ├── assemble_data.py               the derived files of a run, in the layout of data/
│   ├── figstyle.py                    figure style
│   ├── make_figures.py                Figures 1 to 5
│   └── make_supplementary_figures.py  Supplementary Figures S1 and S2
├── units/                             the term lists of the four sub-units, one panel_v2.py each
├── data/                              derived result files: the whole panel at the top, units/ and sensitivity/
└── figures/                           figure1.png to figure5.png, figureS1.png and figureS2.png
```

The ClinVar files, the gnomAD records and the low-complexity intervals are retrieved, not redistributed. `codes/run_all.sh` reads `variant_summary.txt.gz`, `submission_summary.txt.gz` and `gene_condition_source_id` from the repository root and takes the path of `LCR-hs38.bed.gz` as its second argument. `codes/gnomad_pull_durable.py` retrieves the gnomAD v4 records and stores every answer under the work directory, and `codes/fsd2.py` retrieves the RefSeq transcripts it needs into `tx_cache.json`.

## Dataset summary

| File | Rows | Description |
| --- | --- | --- |
| `data/variant_table.csv` | 91,712 | Every pathogenic and likely pathogenic nonsense and frameshift variant of the whole panel's isolated tier and every gnomAD v4 nonsense and frameshift variant of the same genes called high-confidence by LOFTEE and passing the site filters, with gene, class, codon, relative position, reported phenotype side, frameshift-to-stop distance, allele count and escape-window membership. |
| `data/panel_v2.csv` | 1,376 | The genes of the whole panel with tier, ClinVar symbols, UniProt accession, protein length and pathogenic variant counts by class. |
| `data/units/*/variant_table.csv.gz` | 14,700 to 31,718 | The same table for each of the four sub-units. |
| `data/units/*/panel_v2.csv` | 276 to 581 | The genes of each sub-unit. |

The JSON files beside them carry the per-gene results of each test: `derived_frame.json` for the frame tests, `derived_fsd2.json` and `derived_fsd_q.json` for the frameshift-to-stop distance, `derived_pilot_windows.json`, `derived_labels_submission.json` and `derived_position_q.json` for the position axis under the broad rule and the files ending `_strictdual` under the strict rule, and `derived_pooled_exclusions.json` for every pooled test beside the same test with ATM and with the carrying gene removed. `results_summary.json` summarises every unit and restriction run, `gate_brain.json` the forty panel cuts and `grouping_names.json` the grouping names. Under `data/sensitivity/` are the 16 restriction runs and the two grouping alternatives, with the interval file's source and checksum in `lcr_source.json`.

## Reproducing the analysis

Run from the repository root.

```
bash codes/run_all.sh work LCR-hs38.bed.gz
python3 codes/make_figures.py
python3 codes/make_supplementary_figures.py
```

`codes/run_all.sh` writes the gate, the cut ClinVar files, the gnomAD pull and one run directory per unit and restriction under `work/`, and assembles their derived files in `work/data/` in the layout of `data/`. `codes/gnomad_pull_durable.py` and `codes/fsd2.py` are the steps that reach the network, and both resume from the answers already stored. The two figure scripts read `data/`, so a new run is drawn by copying `work/data/` over it. `codes/uniprot_rows.py` rebuilds `codes/uniprot_rows.json` from the UniProt answers of the day it is run.

## Software dependencies

| Domain | Packages |
| --- | --- |
| Statistics | `numpy`, `scipy` |
| Plotting | `matplotlib` |
| Retrieval | Python standard library only |

The analysis ran under Python 3.13.18 with NumPy 2.5.3 and SciPy 1.18.1, and the figures were drawn under Python 3.10.12 with NumPy 2.2.6 and Matplotlib 3.10.9.

## Data availability

* ClinVar: `gene_condition_source_id` dated 9 November 2026, `variant_summary.txt.gz` dated 9 November 2026 and `submission_summary.txt.gz` dated 8 November 2026, from `https://ftp.ncbi.nlm.nih.gov/pub/clinvar/`.
* gnomAD v4 predicted loss-of-function variants: the gnomAD API, `https://gnomad.broadinstitute.org/api`, retrieved 11 November 2026.
* Low-complexity intervals: `scripts/LCR-hs38.bed.gz` of `https://github.com/lh3/varcmp` at commit `bb5b616`.
* UniProt: `https://rest.uniprot.org/uniprotkb/search`, the answers of 11 November 2026 for 1,415 symbols, in `codes/uniprot_rows.json`.
* RefSeq transcripts: NCBI E-utilities, `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi`.

## Citation

Karagöl T., Karagöl A. *Nonsense-mediated decay escape in neurological genes and reported phenotype of truncating variants* (2026).

## Contact

Taner Karagöl, karagolt@mit.edu, ORCID 0009-0005-1011-7661

Alper Karagöl, akaragol@mit.edu, ORCID 0009-0001-7864-0732

## License

The code in this repository is released under the MIT License. The ClinVar, gnomAD, UniProt and RefSeq records it retrieves are subject to the terms of their providers.
