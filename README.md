# Meneval : Meneco validation
Pipeline to run Meneco gapfilling tool and validate reactions proposed

### Meneco Tool :

- Uses Meneco for GapFilling problem solving

S. Prigent et al., “Meneco, a Topology-Based Gap-Filling Tool Applicable to Degraded Genome-Wide Metabolic Networks,” PLOS Computational Biology, vol. 13, no. 1, p. e1005276, Jan. 2017. https://doi.org/10.1371/journal.pcbi.1005276

GitHub : https://github.com/bioasp/meneco

### Padmet Tool :

- Uses padmet networks format to manage all networks and padmets-utils functions.

Aite, M., Chevallier, M., Frioux, C., Trottier, C., Got, J., Cortés, M. P., Mendoza, S. N., Carrier, G., Dameron, O., Guillaudeux, N., Latorre, M., Loira, N., Markov, G. V., Maass, A., and Siegel, A. (2018). Traceability, reproducibility and wiki-exploration for “à-la-carte” reconstructions of genome-scale metabolic models. PLOS Computational Biology, 14(5), e1006146. https://doi.org/10.1371/journal.pcbi.1006146

GitHub : https://github.com/AuReMe/padmet

## Requirements

### Tools
- Must have ncbi-blast installed
  - with `/bin/blastp` 
  - with `/bin/tblastn`

### Python packages

- Pip packages required :
  - `biopython >= 1.80`
  - `Meneco >= 2.0.2`


- Cloned packages required :
  - `aucomana` : https://github.com/PaulineGHG/aucomana.git
  - `padmet` : https://github.com/AuReMe/padmet.git

## Installation

### NCBI-BLAST installation guide 

https://www.ncbi.nlm.nih.gov/books/NBK569861/

### Python packages

From this cloned repository :

```commandline
bash install_dependencies.sh
```

## Usages

**Initialisation of project directories**
```commandline
meneval --init
```
Fill directories with the rights input files

**Checking for input files**
```commandline
meneval --check
```
**Generate supplementary input files**
```commandline
meneval --files
```
**Gapfilling + Adding reactions with blastp hit**
```commandline
meneval --blastp
```
**Gapfilling + Adding reactions from enrichment networks**

```commandline
meneval --enrich Group
```
(`meneval --enrich all` for considering all groups at once)

**Gapfilling + Adding reactions left**
```commandline
meneval --fill
```
**Removing reactions from enrichment networks**
```commandline
meneval --exclude
```
**Run all steps from check step**
```commandline
meneval --workflow
```

## Documentation

See the full documentation on wiki :
https://github.com/PaulineGHG/Meneval/wiki
