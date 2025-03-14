nextflow.enable.dsl=2
// Example command:
// module load python/3.9.3_anaconda2021.11_mamba
// nextflow run main_AMR++.nf -profile conda --pipeline demo
// nextflow run main_AMR++.nf -profile conda --pipeline demo --kraken_db /mnt/c/Users/enriq/Dropbox/minikraken_8GB_20200312/

/*
 * Default pipeline parameters. They can be overriden on the command line eg.
 * given `params.foo` specify on the run command line `--foo some_value`.
*/

log.info """\
 A M R + +    N F   P I P E L I N E
 ===================================
 reads        : ${params.reads}
 pipeline     : ${params.pipeline}
 output       : ${params.output}
 """


def helpMessage = """\
    AMR++ Nextflow Pipeline Help
    =============================

    Available pipelines:
        - demo: Run a demonstration of AMR++
        - standard_AMR: Run the standard AMR++ pipeline
        - fast_AMR: Run the fast AMR++ pipeline without host removal.
        - standard_AMR_wKraken: Run the standard AMR++ pipeline with Kraken
    Available pipeline subworkflows:
        - eval_qc: Run FastQC analysis
        - trim_qc: Run trimming and quality control
        - rm_host: Remove host reads
        - resistome: Perform resistome analysis
        - align: Perform alignment to MEGARes database
        - kraken: Perform Kraken analysis
        - qiime2: Perform QIIME 2 analysis
        - bam_resistome: Perform resistome analysis on BAM files

    To run a specific pipeline/subworkflow, use the "--pipeline" option followed by the pipeline name:
        nextflow run main_AMR++.nf --pipeline <pipeline_name> [other_options]

    To analyze your samples or otherwise change how AMR++ runs, modify the "params.config" file 
    or add more parameters to the command line.

    Finally, consider your computing environment and modify the "-profile" option. By default,
    AMR++ assumes all software dependencies are in your \$PATH, as in the "local" profile. Here are 
    the other options:
        - local: Assumes all sofware is already in your \$PATH
        - local_slurm: Local installation and adds control over slurm job submission.
        - conda: Uses "mamba" to install a conda environment. 
        - conda_slurm: Uses "mamba" and adds control over slurm job submission.
        - singularity: Uses a "singularity" image container.
        - singularity_slurm: Singularity image and adds control over slurm job submission.
        - docker: Uses a docker image container.

    To include SNP analysis, add `--snp Y` to your command.
    
    To include deduplicated count analysis, add `--deduped Y` to your command. 
    Please be aware that adding deduplicated counts will significantly increase run time and temp file storage requirements.

    """

Channel
    .fromFilePairs( params.reads , size: (params.reads =~ /\{/) ? 2 : 1)
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .map { id, files -> 
        def modified_baseName = id.split('\\.')[0]
        tuple(modified_baseName, files)
    }
    .set {fastq_files}

// Load null pipeline
params.pipeline = null

// Load main pipeline workflows

// Load subworkflows
include { FASTQ_QC_WF } from './subworkflows/fastq_information.nf'
include { FASTQ_TRIM_WF } from './subworkflows/fastq_QC_trimming.nf'
include { FASTQ_ALIGN_WF } from './subworkflows/fastq_align.nf'
include { FASTQ_ALIGN_TO_ALL_WF } from './subworkflows/fastq_align_all_to_all.nf'
include { FASTQ_ONLY_ALIGN_TO_ALL_WF } from './subworkflows/fastq_only_align_all_to_all.nf'
include { FASTQ_RM_HOST_WF } from './subworkflows/fastq_host_removal.nf' 
include { FASTQ_KRAKEN_WF } from './subworkflows/fastq_microbiome.nf'
include { FASTQ_SKESA_WF } from './subworkflows/fastq_assembly.nf' // might erase
include { FASTQ_KRAKEN_EXTRACT_WF } from './subworkflows/fastq_microbiome_only_extract.nf'
include { FASTQ_KRAKEN_SINGLE_SPECIES_WF } from './subworkflows/fastq_microbiome_single_extract_species.nf'
include { FASTQ_KRAKEN_DOUBLE_SPECIES_WF } from './subworkflows/fastq_microbiome_double_extract_species.nf'
include { FASTQ_KRAKEN_ONLY_CONFIRMATION_WF } from './subworkflows/fastq_microbiome_only_confirmation.nf'
include { FASTQ_PSEUDOALIGN_WF } from './subworkflows/fastq_pseudoalign.nf'
include { FASTQ_MSWEEP_WF } from './subworkflows/fastq_pseudoalign.nf'
include { FASTQ_DEDUP_WF } from './subworkflows/fastq_dedup.nf'
include { FASTQ_DEDUP_BBMAP_WF } from './subworkflows/fastq_dedup_bbmap.nf'


workflow {
    if (params.pipeline == null || params.pipeline == "help") {

        println helpMessage


        log.info """\
        ===================================
        Running a demonstration of AMR++
        ===================================
        """
        //run with demo params, use params.config
        FASTQ_DEDUP_WF(fastq_files)
        
    }
    else if(params.pipeline == "demo") {
        log.info """\
        ===================================
        Running a demonstration of AMR++
        ===================================
        """
        //run with demo params, use params.config
        FASTQ_DEDUP_WF(fastq_files)
    } 
    else if(params.pipeline == "dedup_cdhit") {

        FASTQ_DEDUP_WF( fastq_files )
    } 
    else if(params.pipeline == "dedup") {

        FASTQ_DEDUP_BBMAP_WF( fastq_files )
    } 
    else if(params.pipeline == "eval_qc") {

        FASTQ_QC_WF( fastq_files )
    } 
    else if(params.pipeline == "trim_qc") {

        FASTQ_TRIM_WF( fastq_files )
    }
    else if(params.pipeline == "rm_host") {

        FASTQ_RM_HOST_WF(params.host, fastq_files )
    } 
    else if(params.pipeline == "align") {

        FASTQ_ALIGN_WF( fastq_files, params.amr)
    }
    else if(params.pipeline == "align_to_all") {

        FASTQ_ALIGN_TO_ALL_WF( fastq_files, params.amr)
    }  
    else if(params.pipeline == "kraken") {
       FASTQ_KRAKEN_WF(fastq_files, params.kraken_db)
    }
    else if(params.pipeline == "only_extract") {
       FASTQ_KRAKEN_EXTRACT_WF(fastq_files, params.kraken_db)
    }
    else if(params.pipeline == "only_confirmation") {
       FASTQ_KRAKEN_ONLY_CONFIRMATION_WF(fastq_files, params.confirmation_db)
    }
    else if(params.pipeline == "single_extract") {
       FASTQ_KRAKEN_SINGLE_SPECIES_WF(fastq_files, params.kraken_db, params.confirmation_db)
    }
    else if(params.pipeline == "double_extract") {
       FASTQ_KRAKEN_DOUBLE_SPECIES_WF(fastq_files, params.kraken_db, params.krakendb_inter, params.confirmation_db)
    }
    else if(params.pipeline == "assembly") {
        FASTQ_SKESA_WF( fastq_files )
    }
    else if(params.pipeline == "only_align_to_all") {
        FASTQ_ONLY_ALIGN_TO_ALL_WF( fastq_files, params.genome_ref_dir)
    } 
    else if(params.pipeline == "pseudoalign") {
        FASTQ_PSEUDOALIGN_WF ( fastq_files)
    }
    else if(params.pipeline == "msweep") {
        FASTQ_MSWEEP_WF( fastq_files )
    }

    else if(params.pipeline == "qiime2") {
        Channel
            .fromFilePairs( params.reads, flat: true )
            .ifEmpty { exit 1, "Read pair files could not be found: ${params.reads}" }
            .map { name, forward, reverse -> [ forward.drop(forward.findLastIndexOf{"/"})[0], forward, reverse ] } //extract file name
            .map { name, forward, reverse -> [ name.toString().take(name.toString().indexOf("_")), forward, reverse ] } //extract sample name
            .map { name, forward, reverse -> [ name +","+ forward + ",forward\n" + name +","+ reverse +",reverse" ] } //prepare basic synthax
            .flatten()
            .collectFile(name: 'manifest.txt', newLine: true, storeDir: "${params.output}/demux", seed: "sample-id,absolute-filepath,direction")
            .set { ch_manifest }
        
        FASTQ_QIIME2_WF( ch_manifest , params.dada2_db)
    }
    else {
            println "ERROR ################################################################"
            println "Please choose a pipeline!!!" 
            println ""
            println "To test the pipeline, use the \"demo\" pipeline or omit the pipeline flag:"
            println ""
            println "ERROR ################################################################"
            println helpMessage
            println "Exiting ..."
            System.exit(0)  
    }
}


workflow.onComplete {
    println "Pipeline completed!"
    println "Started at  $workflow.start"
    println "Finished at $workflow.complete"
    println "Time elapsed: $workflow.duration"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}
