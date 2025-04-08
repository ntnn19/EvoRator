# EvoRator
EvoRator

## Clone this repository
```bash
git clone --recurse-submodules https://github.com/ntnn19/EvoRator.git
```

## Setup the virtual environment 
```bash
cd EvoRator
micromamba env create -p $PWD/venv -f env.yaml
```
 
## Run an example
```bash
python evorator.py 2lzm "" A \
data/catalytic_sites.csv output --orphan-prediction="True" \
 --trained-regressor=evorator_model/EvoRator_SVR_final_model_PERFO_281121.joblib \
--consurfdb-query="" --consurf-output="" --prediction-task="" --job-title=2LZM
```
