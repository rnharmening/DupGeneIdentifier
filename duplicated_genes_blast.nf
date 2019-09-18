#!/usr/bin/env nextflow

/*
 *
 * Based on an example by the nextflow author(s).
 * 
 */
 
/*
 * The pipeline inputs parameters which have to be specified as command line options
 */

blast_cmd = params.blast_cmd
blastFile = file(params.blastFile)
blast_dbName = params.blast_dbName
dbDir = file(params.dbDir)
pubDir = file(params.pubDir)

/*
 * optional input paramters
 */

params.max_target_seqs = 100
params.pident_threshold = 80
params.qcovs_threshold = 75


dbtypes = ["blastn":"nucl", "blastp":"prot"] //, "blastx":"prot", "tblastn":"nucl", "tblastx":"nucl"]

dbtype = dbtypes[blast_cmd]
 

println "Searching for duplicated genes in " + blastFile + " using itself as the database"
println "results stored in: " + pubDir


// make publishDir if it doesn't exist
if( !pubDir.exists() ) {
  new File("$pubDir").mkdir()  
}

/*
 * 1a) Make BLAST database for the input file
 */
process makedb {
    storeDir "$dbDir"

    input:
    file blastFile

    output:
    //file db
    file "${blast_dbName}.*" into db

    """
    makeblastdb -in $blastFile -dbtype $dbtype -out $blast_dbName -parse_seqids -hash_index 
    """
}

/* 
 * 2) blast file against its database
 */
process blast {
    publishDir "$pubDir", mode: 'copy'

    input:
    file db
    file blastFile
    val dbDir
    val blast_dbName

    output:
    file out_file into blast_table_ch

    script:
    out_file = "${blastFile.simpleName}.blast_hits.tsv"
    """
    $blast_cmd -query $blastFile \
    -db $dbDir/$blast_dbName \
    -out $out_file \
    -outfmt '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qframe sframe qlen slen qcovs qcovhsp' \
    -evalue 0.0001 \
    -max_target_seqs $params.max_target_seqs -num_threads 10 
    """
}

/*
 * 3) run python script bin/identify_duplicated_genes.py
 */
process identify_duplicated_genes {
    publishDir "$pubDir", mode: 'copy'

    input:
    file blast_table from blast_table_ch
 
    output:
    file out_file

    script:
    out_file = "${blast_table.simpleName}.duplicates.tsv"
    """
    identify_duplicated_genes.py \
        -i $blast_table \
        -o $out_file \
        --pident $params.pident_threshold \
        --qcovs $params.qcovs_threshold
    """
}
  
workflow.onComplete {
      def subject = 'Duplicated Gene Identification'
      
      def msg = """\
      Pipeline execution summary  
      ---------------------------
      Cmd line    : ${workflow.commandLine}
      Completed at: ${workflow.complete}
      Duration    : ${workflow.duration}
      Success     : ${workflow.success}
      workDir     : ${workflow.workDir}
      publishDir  : $pubDir
      exit status : ${workflow.exitStatus}
      Error report: ${workflow.errorReport ?: '-'}
      """
      .stripIndent()

      println(msg)
  }
