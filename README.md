# biogeographic-modelling-of-arctic-dogs

This repository accompanies Josh Harkness' MSc dissertation project (supervised by Davide Pisani) on Arctic dog domestication using mitogenome data. It aggregates the inputs, analysis scripts, outputs, and metadata required to reproduce the BioGeoBEARS comparative biogeographic analyses, stochastic mapping, and supporting phylogenetic work.

## Abstract
The dispersal of domestic dogs (Canis lupus familiaris) into the Americas has long been  38 debated, particularly whether Arctic lineages reflect a single Thule-era replacement or a  39 series of Holocene introductions. In this study, I combine mitochondrial phylogenomics with  40 biogeographic modelling to reconstruct the colonisation history of Arctic dogs and examine  41 their role in human–dog co-dispersal. The dataset comprised mitochondrial genomes from  42 dogs, wolves, and coyotes, along with newly assembled Greenland sled dog sequences.  43 Phylogenetic inference in IQ-TREE, ancestral range estimation in BioGeoBEARS, and  44 stochastic mapping to estimate dispersal frequencies collectively suggest a complex history  45 rather than a simple introduction. A fossil-constrained, time-scaled phylogeny provided 46 temporal context for these events; placing repeated trans-Beringian dispersals at distinct  47 points during the Holocene. Eastern Siberian populations consistently appear as sources for  48 Arctic lineages, suggesting multiple entry waves. These findings support a layered model of  49 canine dispersal, mirroring human movement, and challenge simplified accounts of total 50 replacements. 

## Repository Structure
- `assemblies/` — consensus mitogenome assemblies and supporting scripts/reference files.
- `alignments/` — per-gene MAFFT alignments, concatenated matrices, and alignment utilities.
- `Figures+Tables/` — publication-ready figures (model comparison, BSM summaries, phylogeny renders).
- `IQTREE/` — phylogenetic inference outputs (`*.treefile`, `.iqtree`, logs) and the `run_iqtree` helper.
- `inputs/` — BioGeoBEARS inputs (e.g. rooted tree, geographic range file).
- `metadata/` — sample-level spreadsheets (dog metadata, SRA list) for provenance tracking.
- `outputs/` — model results (`BioGeoBEARS_Results/`) and BSM artefacts (`BSM_Results/`).
- `R analysis/` — R scripts for BioGeoBEARS model batches, BSM runs, plotting, and time-scaling.

## Software Requirements
- R (>= 4.0 recommended) with packages listed in `R analysis/requirements.txt`.
- BioGeoBEARS and dependencies (`ape`, `phangorn`, `dplyr`, `ggplot2`, `ggtree`, etc.).
- Optional: IQ-TREE 2 for tree inference (`IQTREE/`).

Install R packages via:

```r
pkgs <- scan("R analysis/requirements.txt", what = character())
missing <- setdiff(pkgs, installed.packages()[, "Package"])
if (length(missing)) install.packages(missing)
```

## Reproducing Analyses
1. **Phylogeny (optional rerun)**
   - Input alignments: `alignments/`.
   - Run IQ-TREE with the provided `run_iqtree` script if tree inference needs to be repeated.
2. **BioGeoBEARS models**
   - Primary scripts: `R analysis/BioGeoBEARS/dec_batch.R`, `divalike_batch.R`, `bayarealike_batch.R`.
   - SLURM helpers: `run_DEC_models.sh`, `run_divalike.slurm`, `run_bayarealike.slurm` (BluePebble HPC).
   - Inputs expected in `inputs/` (tree: `FcC_supermatrix.coyote1.rooted.fixed.tre`; geography: `biogeobears_input_2.geog`).
3. **Biogeographical Stochastic Mapping (BSM)**
   - Scripts under `R analysis/BSM/` (`BSM.R`, `BSM2.R`, `BSM3.R`, `run_BSM.R`, `run_BSMBAY.R`).
   - Outputs stored in `outputs/BSM_Results/` (e.g. DECJ 500-replicate results and events tables).
4. **Plotting & Summaries**
   - Visualization scripts in `R analysis/Plotting/` produce ancestral state pies, event summaries, and model comparisons found in `Figures+Tables/`.

## Key Outputs
- Model fit objects for DEC, DEC+J, DIVALIKE, DIVALIKE+J, BAYAREALIKE, BAYAREALIKE+J (`outputs/BioGeoBEARS_Results/`).
- BSM replicates and event tables (`outputs/BSM_Results/DECJ_*`, `BAYAREALIKEJ_*`).
- Final figures summarising dispersal scenarios (`Figures+Tables/`).

## Metadata & Provenance
- Sample metadata: `metadata/dog_metadata_filtered_unique.xlsx` and `metadata/sra_list.txt`.
- HPC runtime logs: `outputs/BSM_Results/BSM.out`, `outputs/BSM_Results/BSMBAY.out`.
- Additional provenance artefacts (e.g. session info, SLURM stderr) should be added alongside the relevant outputs.

## Documentation Guidance
This top-level README gives the project overview. If future collaborators need folder-specific instructions (e.g. detailed alignment workflows or figure generation steps), add concise `README.md` files within the corresponding subdirectories. A good rule: create a per-folder README when the scripts or data inside require more than a short paragraph to explain or when reproducing them involves multiple steps or toolchains.

## Citation
If you use this repository, please cite the forthcoming publication and acknowledge BioGeoBEARS, IQ-TREE, and other tools as appropriate.
