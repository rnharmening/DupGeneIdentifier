# DupGeneIdentifier
Nextflow script to identify duplicated gene sequences in a FastA file. Each gene has to be one sequence in the Fasta file

To run the pipeline blast has to be installed.

General Usage:

```bash
nextflow run duplicated_genes_blast.nf \
  --blast_cmd [blastp, blastn] \
  --blastFile <InputFastaFile>\
  --pubDir <ResultsDir> \
  --dbDir <blastDBDir> \
  --blast_dbName <DBbaseName>
```
