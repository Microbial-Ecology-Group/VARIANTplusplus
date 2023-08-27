params.taxlevel = "S" //level to estimate abundance at [options: D,P,C,O,F,G,S] (default: S)
params.readlen = 150

threads = params.threads

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
            else if(filename.indexOf(".kraken.filtered.report") > 0) "Kraken/filtered_report/$filename"
            else if(filename.indexOf(".kraken.filtered.raw") > 0) "Kraken/filtered/$filename"
            else {}
        }

    input:
       tuple val(sample_id), path(reads)
       path(krakendb)


   output:
      tuple val(sample_id), path("${sample_id}.kraken.raw"), emit: kraken_raw
      path("${sample_id}.kraken.report"), emit: kraken_report
      tuple val(sample_id), path("${sample_id}.kraken.filtered.raw"), emit: kraken_filter_raw
      path("${sample_id}.kraken.filtered.report"), emit: kraken_filter_report
      tuple val(sample_id), path("${sample_id}_kraken2.krona"), emit: krakenkrona
      tuple val(sample_id), path("${sample_id}_kraken2_filtered.krona"), emit: krakenkrona_filtered



     """
     ${KRAKEN2} --db ${krakendb} --paired ${reads[0]} ${reads[1]} --threads ${threads} --report ${sample_id}.kraken.report > ${sample_id}.kraken.raw
     ${KRAKEN2} --db ${krakendb} --confidence 1 --paired ${reads[0]} ${reads[1]} --threads ${threads} --report ${sample_id}.kraken.filtered.report > ${sample_id}.kraken.filtered.raw

    cut -f 2,3  ${sample_id}.kraken.raw > ${sample_id}_kraken2.krona
    cut -f 2,3  ${sample_id}.kraken.filtered.raw > ${sample_id}_kraken2_filtered.krona
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
            else if(filename.indexOf(".confirmation.kraken.filtered.report") > 0) "Kraken/confirmation_filtered_report/$filename"
            else if(filename.indexOf(".confirmation.kraken.filtered.raw") > 0) "Kraken/confirmation_filtered/$filename"
            else {}
        }

    input:
        tuple val(sample_id), path(extracted_r1), path(extracted_r2)
        path(confirmation_db)

    output:
        tuple val(sample_id), path("${sample_id}.confirmation.kraken.raw"), emit: confirmation_kraken_raw
        path("${sample_id}.confirmation.kraken.report"), emit: confirmation_kraken_report
        tuple val(sample_id), path("${sample_id}.confirmation.kraken.filtered.raw"), emit: confirmation_filtered_kraken_raw
        path("${sample_id}.confirmation.kraken.filtered.report"), emit: confirmation_filtered_kraken_report

    script:
    """
    ${KRAKEN2} --db ${confirmation_db} --paired ${extracted_r1} ${extracted_r2} --threads ${threads} \
        --report ${sample_id}.confirmation.kraken.report > ${sample_id}.confirmation.kraken.raw

    ${KRAKEN2} --db ${krakendb} --confidence 1 --paired ${extracted_r1} ${extracted_r2} --threads ${threads} \ 
        --report ${sample_id}.confirmation.kraken.filtered.report > ${sample_id}.confirmation.kraken.filtered.raw


    """
}


process extractKrakenReads {
    tag { sample_id }
    label "microbiome"

    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3

    publishDir "${params.output}/MicrobiomeAnalysis/ExtractedReads", mode: 'copy'

    input:
        tuple val(sample_id), path(reads)
        path(kraken_raw)
        path(kraken_report)

    output:
        tuple val(sample_id), path("extracted-r1.fastq"), path("extracted-r2.fastq")

    script:
    """
    python extract_kraken_reads.py -k ${kraken_raw} --report ${kraken_report} --taxid 75984 --include-children --include-parents \
        --fastq-output -s1 ${reads[0]} -s2 ${reads[1]} -o extracted-r1.fastq -o2 extracted-r2.fastq
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
        path(kraken_filtered_reports)

    output:
        path("kraken_analytic_matrix.csv")

    """
    ${PYTHON3} $baseDir/bin/kraken2_long_to_wide.py -i ${kraken_reports} -o kraken_analytic_matrix.csv

    ${PYTHON3} $baseDir/bin/kraken2_long_to_wide.py -i ${kraken_filtered_reports} -o kraken_filtered_analytic_matrix.csv
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
        path(confirmation_filtered_kraken_reports)

    output:
        path("Extracted_species_kraken_analytic_matrix.csv")

    script:
    """
    ${PYTHON3} $baseDir/bin/kraken2_long_to_wide.py -i ${confirmation_kraken_reports} -o Extracted_species_kraken_analytic_matrix.csv
    
    ${PYTHON3} $baseDir/bin/kraken2_long_to_wide.py -i ${confirmation_filtered_kraken_reports} -o Extracted_species_kraken_filtered_analytic_matrix.csv

    """
}




process runbracken {
    label "microbiome"
    
    input:
       tuple val(sample_id), path(krakenout)
       tuple val(sample_id), path(krakenout_filtered)
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