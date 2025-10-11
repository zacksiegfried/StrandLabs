<img width="551" height="141" alt="StrandLabs" src="https://github.com/user-attachments/assets/0847f1b5-090a-40ba-8dff-20a308e9bb7a" />

# StrandLabs is a place to analyze DNA

#### Basic set up guide
- Clone repo onto own machine
- Set up virtual environment
```
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```
- in pipeline folder create "data", "output" and "results" sub folders

### Running pipeline with Nextflow
```
nextflow run pipeline/main.nf -profile test
```
