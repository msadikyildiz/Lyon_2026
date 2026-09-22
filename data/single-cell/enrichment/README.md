# EcoCyc enrichment results

Five EcoCyc SmartTable text exports retain their original bytes. `manifest.json` maps four to the source workbook: 188 cluster-4 results, 163 cluster-9 results, 376 results for genes higher in both evolved cultures, and parent-up results. Cluster-1 and cluster-2 enrichment results are workbook inputs without matching text exports. The separate 74-row “4-and-7-combined” export represents an additional comparison.

Supplementary Table 4 and Source Data contain all 133 parent-up export rows with P < 0.1. The source workbook contains the 19-row P < 0.001 subset. Every overlapping term, P-value and matched-gene list agrees between the text exports and workbook. The original Excel source is unchanged. The reporting cutoff of 0.1 includes results above 0.05.

## Method and assumptions

The documented standard workflow uses EcoCyc’s E. coli K-12 MG1655 database, Fisher’s exact test, Benjamini–Hochberg correction and a reporting cutoff of P < 0.1, covering pathways, transcriptional/translational regulators and Gene Ontology terms. The settings screenshot documents this standard workflow rather than saved settings for each result table. Its 319-gene SmartTable was created on 7 May 2025.

[EcoCyc documentation, section 7.6.1](https://ecocyc.org/PToolsWebsiteHowto.shtml) describes over-representation testing. Default reference sets depend on the annotation class: pathway-annotated genes for pathway tests and GO-annotated genes for GO tests. These default backgrounds are assumed for this analysis. Export object labels such as `All-Genes` and `DNA-Segments` identify result objects, not statistical backgrounds. All exported rows and empty Matches fields are preserved.

[Release notes](https://biocyc.org/ecoli/release-notes.shtml) date EcoCyc 29.0 to 10 April 2025 and 29.1 to 22 September 2025. Version 29.0 is assumed from the May 2025 workbook and SmartTable dates. The database version and background assignments are assumptions, not recovered run metadata. Result tables are repository inputs; the reproducibility workflow verifies and exports them without rerunning enrichment tests.
