# Unified pipeline to ingest cell type specific differential expression into the InTraC knowledge base

The pipeline starts with annotated single cell data in anndata objects.
Depending on the number of samples present in the file the pipeline performs either pseudo bulk differential gene expression (DEG) analysis using edgeR (for more than 3 samples per condition) or using the Wilcoxon rank rum test (if less than 3 samples per condition are present).

## Installation
The pipeline depends on the following packages:
- python
- scanpy
- R
- edgeR

## Input
As input it requires meta data information about 
- the anndata file:
  - name of the cell type column in `.obs`
  - name of the sample id column in `.obs`
  - name of the columns that holds the condition information in `.obs`
  - value of the 'reference' condition to compare to
- the data set:
  - species
  - tissue
  - reference to publication
  - publication year

## Execution

```
python pipeline_main.py [args]
```

## Output
The output is a csv file in canonical form for ingestion with the InTraC knowledge graph adapter.
