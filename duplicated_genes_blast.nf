#!/usr/bin/env nextflow

/*
 *
 * Based on an example by the nextflow author(s).
 * 
 */

def helpMessage() {
    log.info"""
    =========================================
    Dupplicated Gene Finder
    =========================================
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run duplicated_genes_blast.nf --blast_cmd <blastp> --blastFile <.faa> --blast_dbName <db>
    
    Mandatory arguments:
      --blast_cmd             Which blast algorithm to use (blastp or blastn)
      --blastFile             Input File in fasta format, containing protein or nucleotide sequences. Each entry (">") is considered as one gene.
      --blast_dbName <db>     Name how the generated database should be stored
    
    Optional arguments:
      --max_target_seqs       The maximal number of hits for blast to report (take care if you use an < 2.9 version of blast!) [def: 100]
      --pident_threshold      The pident threshold for which genes should be reported as dupplicated [def: blastn:90 blastp:80] 
      --qcovs_threshold       The qcovs threshold for which genes should be reported as dupplicated [def: blastn:90 blastp:80] 
      --eval                  The evalue cutoff for blast [def: blastn:1e-7 blastp:1e-4] 
        
      --dbDir                 The directory where the generated database will be stored [def: DataBases]
      --pubDir                The directory where the results will be stored [def: Results]
    """.stripIndent()
}
 
// Show help message
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

/*
 * The pipeline inputs parameters which have to be specified as command line options
 */

blast_cmd = params.blast_cmd
blastFile = file(params.blastFile)
blast_dbName = params.blast_dbName
params.dbDir = "DataBases"
dbDir = file(params.dbDir)
params.pubDir = "Results"
pubDir = file(params.pubDir)

/*
 * optional input paramters
 */

params.max_target_seqs = 100 
params.pident_threshold = (blast_cmd == "blastn") ? 90 : 80 
params.qcovs_threshold = (blast_cmd == "blastn") ? 90 : 80 
params.eval = (blast_cmd == "blastn") ? 10**(-7): 10**(-4)

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
    -evalue $params.eval \
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
