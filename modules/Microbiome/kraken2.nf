params.taxlevel = "S" //level to estimate abundance at [options: D,P,C,O,F,G,S] (default: S)
params.readlen = 150

threads = params.threads

kraken_confidence = params.kraken_confidence

extract_reads_taxid = params.extract_reads_taxid
extract_reads_options_single = params.extract_reads_options_single
extract_reads_options_double = params.extract_reads_options_double


krakendb_inter = params.krakendb_inter
confirmation_db = params.confirmation_db

process dlkraken {
    tag { }
    label "python"

    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3

    publishDir "$baseDir/data/kraken_db/", mode: 'copy'

    output:
        path("minikraken_8GB_20200312/")

    """
        wget ftp://ftp.ccb.jhu.edu/pub/data/kraken2_dbs/minikraken_8GB_202003.tgz
        tar -xvzf minikraken_8GB_202003.tgz

    """
}


process runkraken {
    tag { sample_id }
    label "microbiome"

    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3

    publishDir "${params.output}/MicrobiomeAnalysis", mode: 'copy',
        saveAs: { filename ->
            if(filename.indexOf(".kraken.raw") > 0) "Kraken/standard/$filename"
            else if(filename.indexOf(".kraken.report") > 0) "Kraken/standard_report/$filename"
            else {}
        }

    input:
       tuple val(sample_id), path(reads)
       path(krakendb)


   output:
      tuple val(sample_id), path("${sample_id}.kraken.raw"), emit: kraken_raw
      path("${sample_id}.kraken.report"), emit: kraken_report


     """
     ${KRAKEN2} --db ${krakendb} --confidence ${kraken_confidence} --paired ${reads[0]} ${reads[1]} --threads ${threads} --report ${sample_id}.kraken.report > ${sample_id}.kraken.raw

    """
}

process runkraken_double_extract {
    tag { sample_id }
    label "microbiome"

    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3

    publishDir "${params.output}/MicrobiomeAnalysis", mode: 'copy',
        saveAs: { filename ->
            if(filename.indexOf(".nt.kraken.raw") > 0) "Kraken/nt_raw/$filename"
            else if(filename.indexOf(".nt.kraken.report") > 0) "Kraken/nt_report/$filename"
            else if(filename.indexOf(".family.kraken.report") > 0) "Kraken/family_report/$filename"
            else if(filename.indexOf(".family.kraken.raw") > 0) "Kraken/family_raw/$filename"
            else if(filename.indexOf(".fastq") > 0) "Kraken/double_extracted_reads/$filename"
            else {}
        }

    input:
       tuple val(sample_id), path(reads)
       path(krakendb)
       path(krakendb_inter)


   output:
      tuple val(sample_id), path("${sample_id}.kraken.raw"), emit: kraken_raw
      path("${sample_id}.kraken.report"), emit: kraken_report
      tuple val(sample_id), path("${sample_id}_Mh_extracted_R*.fastq.gz"), emit: extracted_reads

     """
     ${KRAKEN2} --db ${krakendb} --confidence ${kraken_confidence} --paired ${reads[0]} ${reads[1]} --threads ${threads} --report ${sample_id}.nt.kraken.report > ${sample_id}.nt.kraken.raw

    extract_kraken_reads.py -k ${sample_id}.nt.kraken.raw --report ${sample_id}.nt.kraken.report --taxid ${extract_reads_taxid} --include-children --include-parents --fastq-output -s1 ${reads[0]} -s2 ${reads[1]} -o temp_extracted-r1.fastq -o2 temp_extracted-r2.fastq

     ${KRAKEN2} --db ${krakendb_inter} --confidence ${kraken_confidence} --paired temp_extracted-r1.fastq temp_extracted-r2.fastq --threads ${threads} --report ${sample_id}.family.kraken.report > ${sample_id}.family.kraken.raw

    extract_kraken_reads.py -k ${sample_id}.kraken.raw --report ${sample_id}.kraken.report --taxid ${extract_reads_taxid} ${extract_reads_options_double} --fastq-output -s1 temp_extracted-r1.fastq -s2 temp_extracted-r2.fastq -o ${sample_id}_Mh_extracted_R1.fastq -o2 ${sample_id}_Mh_extracted_R2.fastq
    gzip *Mh_extracted*

    rm temp*

    """
}

process runkraken_extract {
    tag { sample_id }
    label "microbiome"

    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3

    publishDir "${params.output}/MicrobiomeAnalysis", mode: 'copy',
        saveAs: { filename ->
            if(filename.indexOf(".kraken.raw") > 0) "Kraken/standard/$filename"
            else if(filename.indexOf(".kraken.report") > 0) "Kraken/standard_report/$filename"
            else if(filename.indexOf(".fastq.gz") > 0) "Kraken/extracted_reads/$filename"
            else {}
        }

    input:
       tuple val(sample_id), path(reads)
       path(krakendb)


   output:
      tuple val(sample_id), path("${sample_id}.kraken.raw"), emit: kraken_raw
      path("${sample_id}.kraken.report"), emit: kraken_report
      tuple val(sample_id), path("${sample_id}_Mh_extracted_R*.fastq.gz"), emit: extracted_reads

     """
    ${KRAKEN2} --db ${krakendb} --confidence ${kraken_confidence} --paired ${reads[0]} ${reads[1]} --threads ${threads} --report ${sample_id}.kraken.report > ${sample_id}.kraken.raw

    extract_kraken_reads.py -k ${sample_id}.kraken.raw --report ${sample_id}.kraken.report --taxid ${extract_reads_taxid} ${extract_reads_options_single} --fastq-output -s1 ${reads[0]} -s2 ${reads[1]} -o ${sample_id}_Mh_extracted_R1.fastq -o2 ${sample_id}_Mh_extracted_R2.fastq

    gzip *Mh_extracted*

    """
}

process runConfirmationKraken {
    tag { sample_id }
    label "microbiome"

    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3

    publishDir "${params.output}/MicrobiomeAnalysis/Confirmation", mode: 'copy',
        saveAs: { filename ->
            if(filename.indexOf(".confirmation.kraken.raw") > 0) "Kraken/confirmation/$filename"
            else if(filename.indexOf(".confirmation.kraken.report") > 0) "Kraken/confirmation_report/$filename"
            else {}
        }

    input:
        tuple val(sample_id), path(extracted_reads)
        path(confirmation_db)

    output:
        tuple val(sample_id), path("${sample_id}.confirmation.kraken.raw"), emit: confirmation_kraken_raw
        path("${sample_id}.confirmation.kraken.report"), emit: confirmation_kraken_report

    script:
    """
    ${KRAKEN2} --db ${confirmation_db} --confidence ${kraken_confidence} --paired ${extracted_reads[0]} ${extracted_reads[1]} --threads ${threads} --report ${sample_id}.confirmation.kraken.report > ${sample_id}.confirmation.kraken.raw

    """
}





process krakenresults {
    tag { }
    label "python"

    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3

    publishDir "${params.output}/Results/", mode: 'copy'

    input:
        path(kraken_reports)

    output:
        path("kraken_analytic_matrix.csv")

    """
    ${PYTHON3} $baseDir/bin/kraken2_long_to_wide.py -i ${kraken_reports} -o kraken_analytic_matrix.csv
    """
}

process extractedKrakenResults {
    tag { }
    label "python"

    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3

    publishDir "${params.output}/Results/", mode: 'copy'

    input:
        path(confirmation_kraken_reports)

    output:
        path("strain_classification_kraken_matrix.csv")
        path("taxonomy_strain_classification_kraken_matrix.csv")
        path("summary_strain_classification_kraken_matrix.csv")
        path("dataset_summary_strain_classification_kraken_matrix.csv")
        

    script:
    """
    ${PYTHON3} $baseDir/bin/kraken2_strains_long_to_wide_wsummary.py -i ${confirmation_kraken_reports} -o strain_classification_kraken_matrix.csv
    

    """
}




process runbracken {
    label "microbiome"
    
    input:
       tuple val(sample_id), path(krakenout)
       path(krakendb)

    """
    bracken \
        -d ${krakendb} \
        -r ${params.readlen} \
        -i ${krakenout} \
        -l ${params.taxlevel} \
        -o ${sample_id}_bracken.tsv

    bracken \
        -d ${krakendb} \
        -r ${params.readlen}\
        -i ${krakenout_filtered} \
        -l ${params.taxlevel} \
        -o ${sample_id}_bracken_filtered.tsv
        """
}

process kronadb {
    label "microbiome"
    output:
        file("krona_db/taxonomy.tab") optional true into krona_db_ch // is this a value ch?

    when: 
        !params.skip_krona
        
    script:
    """
    ktUpdateTaxonomy.sh krona_db
    """
}

process kronafromkraken {
    publishDir params.outdir, mode: 'copy'
    label "microbiome"
    input:
        file(x) from kraken2krona_ch.collect()
        //file(y) from kaiju2krona_ch.collect()
        file("krona_db/taxonomy.tab") from krona_db_ch
    
    output:
        file("*_taxonomy_krona.html")

    when:
        !params.skip_krona
    
    script:
    """
    mkdir krona
    ktImportTaxonomy -o kraken2_taxonomy_krona.html -tax krona_db $x
    """
}