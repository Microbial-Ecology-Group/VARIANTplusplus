
threads = params.threads



process run_skesa_reads {
    tag "$pair_id"
    label "assembly"

    memory { 2.GB * task.attempt }
    time { 1.hour * task.attempt }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3

    publishDir "${params.output}/Assembly/", mode: "copy",
        saveAs: { filename ->
            if(filename.indexOf(".skesa.fa") > 0) "skesa_contigs/$filename"
            else {}
        }

    input:
        tuple val(pair_id), path(reads) 

    output:
        tuple val(pair_id), path("${pair_id}.skesa.fa"), emit: skesa_contigs

    """
    skesa --reads ${reads[0]},${reads[1]} --cores ${threads} --memory ${mem_skesa} > ${pair_id}.skesa.fa

    """
}
