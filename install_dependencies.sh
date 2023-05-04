#!/bin/bash

pip install -r requirements.txt

mkdir local_packages
cd local_packages/
git clone https://github.com/PaulineGHG/aucomana.git
git clone https://github.com/AuReMe/padmet.git

cd aucomana/
pip install -e .
cd ../padmet/
pip install -e .

cd ../../