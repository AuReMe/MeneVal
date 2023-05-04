# Meneco_validation
Pipeline to run Meneco gapfilling tool and validate reactions proposed

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
