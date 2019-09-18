# DupGeneIdentifier

Nextflow script to identify duplicated gene sequences in a FastA file. Each gene has to be one sequence in the Fasta file

To run the pipeline blast has to be installed.

## General Usage

```bash
nextflow run duplicated_genes_blast.nf \
  --blast_cmd [blastp, blastn] \
  --blastFile <InputFastaFile>\
  --pubDir <ResultsDir> \
  --dbDir <blastDBDir> \
  --blast_dbName <DBbaseName>
  [-profile conda]
```

If you have conda installed, you can use the conda profile to automatically download the needed software listed in the `environment.yml` file.

## Demo

In the `demo` directory is one `demo.faa` file containing some aminoacid sequences, where some sequences are duplicates of others (with some minor changes). You can run the pipeline with the demo files to verify that everithing works. The output files should match the files in `./demo/demo_results/`. 

```
nextflow run duplicated_genes_blast.nf --blast_cmd blastp --blastFile demo/demo.faa --pubDir demo_results --dbDir demo_DB --blast_dbName demo
```
