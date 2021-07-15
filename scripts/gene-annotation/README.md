# Gene Annotation Scripts

Lots going on here, but there are two main groups of scripts:

## Annotation of FiNZ-ZNF genes
### Original predictions
1. `generate_profiles.sh`
2. `predict_finz.sh`
3. ``

### Retraining augustus for improved D. rerio annotations
1. `retrain_augustus.sh`

## Quality control:
### Read coverage over predicted genes.
1. `get_busco_beds.py`
2. `extract_flanked_genes.sh`
3. `map_danio_reads.sh`
4. `blast_busco_exons.sh`

### Comparison with Ensembl and RefSeq annotations
1. `overlap_finz.sh`


## Accessory scripts used within other progs:
1. `offset_gffs.py`
2. `unmask.py`

