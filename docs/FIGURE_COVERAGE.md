# Figure coverage

“Rebuilt” means recalculated or plotted from included numerical inputs by `reproduce.py`. “Assembled” means the complete figure is composed from regenerated panels plus the listed static panels. “Reference code” means the notebook is included but its complete input data are missing. The detailed [panel manifest](../working/analysis/stats-rework/out/panel_manifest.csv) links each measure to a source; it describes experimental records, while this table states executable coverage.

| Figure | Rebuilt by the command | Retained assets or remaining dependency |
|---|---|---|
| 1 | Survival d–f panels and Source Data | Schematics a–c supplied as PPTX; complete figure composition remains manual |
| 2 | IC50 b–d and MDK e; complete assembly | Representative curve a and genomic f–h retained as images; curve data/cache and genomic reference notebooks plus complete displayed endpoint tables included |
| 3 | IC50 b–d and MDK e; complete assembly | Survival a retained as image in the assembly; its original plotting notebook also executes separately |
| 4 | No figure rendering | Genomic a–c reference notebook needs full PLAC calls/frequency tables; schematic d supplied as PPTX and PDF |
| 5 | Growth a, IC50 b–d, MDK e; complete assembly | Biological n=1 per strain; technical series retained |
| Supplementary 1 | No figure rendering | Original schematic PPTX supplied |
| Supplementary 2 | IC50 a and MIC b; complete assembly | None |
| Supplementary 3 | IC50 b–d; complete assembly | Resistance a and genomic e retained as images in the assembly; resistance plotting notebook also executes separately; complete displayed genomic table supplied |
| Supplementary 4 | IC50 b–d and MDK e; complete assembly | Survival a and genomic f retained as images in the assembly; survival plotting notebook also executes separately; genomic table has 50/83 rows |
| Supplementary 5 | Tables exported to Source Data | Genomic reference notebooks plus complete displayed tables; full upstream calls still needed |
| Supplementary 6 | IC50 a–c; complete assembly | Genomic d retained as image; complete displayed table and reference notebook supplied |
| Supplementary 7 | IC50 a and MIC b; complete assembly | None |
| Supplementary 8 | No figure rendering | Reference notebook needs full PLAC traced-allele data |
| Supplementary 9 | No figure rendering | Reference notebook needs full mutation tables |
| Supplementary 10 | MIC a–c, MDK d; complete assembly | None |
| Supplementary 11 | Recalculated fractions and figure | Absolute fractions conditional on equal plating volumes |
| Supplementary 12–13 and Figure 6 | Source Data dependency entries | Single-cell source objects, metadata, full tables and code missing |
| Supplementary Tables 1–3 | Source Data sheets rebuilt from supplied JSON | Original table values extracted from the manuscript |
| Supplementary Table 4 | Source Data dependency entry | Full single-cell source table and mapping needed |

The ten assemblies contain 50 placements: 40 regenerated size-specific panels and ten retained images. `assembly_assets.json` maps all input assets and `assembly_layout.json` records their layout; `out/print-panels/manifest.json` identifies replacements generated from numerical inputs. The output assembly manifest records final placements and validation.

The six original trajectory notebooks also render standalone panels to `out/trajectories`. These retain the original plotting rules, with the validated DejaVu Sans fallback font selected explicitly in the survival notebooks; the assembled figures continue to use the established static trajectory panels.
