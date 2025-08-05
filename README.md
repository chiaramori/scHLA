# scHLA
Analysis of sequencing results from scHLA/phage display experiments with the goal of identifying what peptides bind the target TCR

There are three R files but only two need to be actuually run:

## -> scHLA_01_WebLogo&Motif
### Description:
- Processes 11aa insert sequences from barcode folders
- Trims and saves to `trimmed/`
- Plots sequence logos
- Classifies by AA class motif
- Computes motif frequencies and plots barplots
- Computes and plots dendrograms of motifs

## -> scHLA_02_ModifyDatabase
- Clean up for database downloaded from IEDB (epitope_full_v3.tsv)
- Changed headings and selected only relevant columns (aa sequence, source_protein and source_organism)
- If peptide belongs to an analog then source_protein and source_organism are taken from the analog column

## -> scHLA_03_tSNE
- First compare peptides contained in trimmed folder (multiple files allowed) with database created with Script2
- Annotate your files with source_protein and source_organism and save in trimmed_blasted folder
- Perfomr distance computation for clustering and save in distance folder
- Perfomr tsne analysis and save in tsne folder.
- Create tsne plots based on counts (highest count forward), on species and on human/non-human. save in plots folder
- Perfomr umap analysis and save in umap folder.
- Create umpa plots based on counts (highest count forward), on species and on human/non-human. save in plots folder
