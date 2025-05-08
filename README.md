![screenshot](assets/RainStreamLogo2.png)

# RainStream is a program to statistically analyze DNA variant files

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

## p-Score Algorithm
