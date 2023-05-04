# Meneval : Meneco validation
Pipeline to run Meneco gapfilling tool and validate reactions proposed

### Meneco Tool :
S. Prigent et al., “Meneco, a Topology-Based Gap-Filling Tool Applicable to Degraded Genome-Wide Metabolic Networks,” PLOS Computational Biology, vol. 13, no. 1, p. e1005276, Jan. 2017. https://doi.org/10.1371/journal.pcbi.1005276

GitHub : https://github.com/bioasp/meneco

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

```commandline
meneval --init
```
```commandline
meneval --check
```
```commandline
meneval --files
```
```commandline
meneval --blastp
```
```commandline
meneval --holobiont
```
```commandline
meneval --aucome --group=group
```
```commandline
meneval --fill
```
```commandline
meneval --workflow
```

## Documentation

See the full documentation on wiki :
https://github.com/PaulineGHG/Meneval/wiki
