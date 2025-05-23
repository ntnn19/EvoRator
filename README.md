# EvoRator
EvoRator

## Clone this repository
```bash
git clone https://github.com/ntnn19/EvoRator.git
```

## Install dssp
```bash
sudo apt-get install dssp
```

## Setup & activate the virtual environment 
```bash
cd EvoRator
micromamba env create -p $PWD/venv -f env.yaml
micromamba activate $PWD/venv
```
 
## Run an example with a pdb identifier
```bash
python evorator.py 2lzm "" A \
data/catalytic_sites.csv output --orphan-prediction="True" \
 --trained-regressor=evorator_model/EvoRator_SVR_final_model_PERFO_281121.joblib \
--consurfdb-query="" --consurf-output="" --prediction-task="" --job-title=2LZM
```

## Run an example with a pdb file
```bash
python evorator.py "" example/alphafold.pdb \
A data/catalytic_sites.csv output --orphan-prediction="True" \
--trained-regressor=evorator_model/EvoRator_SVR_final_model_PERFO_281121.joblib \
--consurfdb-query="" --consurf-output="" --prediction-task="" --job-title=alphafold
```
