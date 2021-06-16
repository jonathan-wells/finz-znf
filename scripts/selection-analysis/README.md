# Selection Analyses

These scripts are primarily for the purpose of extracting clades of paralogous 
FiNZ-ZnF genes that are tested for signatures of positive selection. The
expression of these clades (in D. rerio) is also checked, and are used to assist
in the design of guide-RNAs for use in dCas9 CRISPRi knockdown experiments.

## Workflow
1. `calc_pairwise_matrix.sh`
    This script calculates all pairwise alignments of FiNZ-ZnF proteins
    sequences from predicted gene annotations, using the Needleman-Wunsch
    algorithm for global alignment. This can easily be adapted to do local
    alignments instead (Smith-Waterman algorithm), but this is likely to yield
    many more spurious relationships due to the highly repetitive nature of the
    proteins. Only those species that pass the quality control checks in
    `analysis/genome_quality_control.ipynb` are used.
2. `parse_pairwise_alignments.py`
    After alignments have been generated, we then calculate the distance between
    sequences using the Jukes-Cantor correction (for proteins) and output the
    resulting data in a Phylip formatted distance matrix.
3. `calc_njtrees.sh`
    Calculates the neighbour-joining tree from the distance matrix just
    mentioned and outputs tree in Newick format.
4. `njtree_labels.py`
    Simple script that replaces leaf labels with more informative versions.
5. `extract_clades.py`
    For the purposes of detecting selection it is important that our genes are
    strictly paralogous or orthologous (as an aside, this may not be such an
    important assumption in practice). We therefore search the resulting NJ tree
    for clades of proteins that are derived from a single species, and for
    practical purposes, contain at least 10 members that differ by no more than
    30 residues in length.
6. `codeml.ctl`
    This is a template for the PAML parameter file used for hypothesis tests of
    selection acting on each clade. Somewhat confusingly, the specific files are 
    kept in the data directory `data/selection-analysis/<clade>`.
