nextflow.enable.dsl=2
// Example command:
// module load python/3.9.3_anaconda2021.11_mamba
// nextflow run main_VARIANT++.nf -profile conda --pipeline demo
// nextflow run main_VARIANT++.nf -profile conda --pipeline demo --kraken_db /mnt/c/Users/enriq/Dropbox/minikraken_8GB_20200312/

/*
 * Default pipeline parameters. They can be overriden on the command line eg.
 * given `params.foo` specify on the run command line `--foo some_value`.
*/

log.info """\
 VARIANT + +    N F   P I P E L I N E
 ===================================
 pipeline     : ${params.pipeline}
 output       : ${params.output}
 """


def helpMessage = """\
    VARIANT++ Nextflow Pipeline Help
    =============================

    Available pipelines:
        - demo: Run a demonstration of VARIANT++
        - full_GSV_pipeline: Run the standard VARIANT++ pipeline
    Available pipeline subworkflows:
        - GSV_1: QC trimming, and merging of reads (--reads)
        - GSV_2: deduplication of reads (--merged_reads)
        - GSV_3: Host DNA removal (--merged_reads)
        - GSV_4: Read filtration with Kraken2 (--merged_reads)
        - GSV_5: GSV Classification with themisto/mSweep (--merged_reads)

    To run a specific pipeline/subworkflow, use the "--pipeline" option followed by the pipeline name:
        nextflow run main_VARIANT++.nf --pipeline <pipeline_name> [other_options]

    To analyze your samples or otherwise change how VARIANT++ runs, modify the "params.config" file 
    or add more parameters to the command line.

    Finally, consider your computing environment and modify the "-profile" option. By default,
    VARIANT++ assumes all software dependencies are in your \$PATH, as in the "local" profile. Here are 
    the other options:
        - local: Assumes all sofware is already in your \$PATH
        - local_slurm: Local installation and adds control over slurm job submission.
        - conda: Uses "mamba" to install a conda environment. 
        - conda_slurm: Uses "mamba" and adds control over slurm job submission.
        - singularity: Uses a "singularity" image container.
        - singularity_slurm: Singularity image and adds control over slurm job submission.
        - docker: Uses a docker image container.


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
include { FLASH_MERGE_WF } from './subworkflows/GSV_flash_reads.nf'
include { GSV_PIPELINE_WF } from './subworkflows/GSV_full_pipeline.nf'
include { GSV_1_WF } from './subworkflows/GSV_step_1_qc_merge.nf'
include { GSV_2_WF } from './subworkflows/GSV_step_2_dedup.nf'
include { GSV_3_WF } from './subworkflows/GSV_step_3_host_rm.nf'
include { GSV_4_WF } from './subworkflows/GSV_step_4_kraken_extraction.nf'
include { GSV_5_WF } from './subworkflows/GSV_step_5_GSV_classification.nf'



workflow {
    if (params.pipeline == null || params.pipeline == "help") {

        println helpMessage


        log.info """\
        ===================================
        Running a demonstration of VARIANT++
        ===================================
        """
        //run with demo params, use params.config
        FASTQ_DEDUP_WF(fastq_files)
        
    }
    else if(params.pipeline == "demo") {
        log.info """\
        ===================================
        Running a demonstration of VARIANT++
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
    else if(params.pipeline == "full_GSV_pipeline") {
        GSV_PIPELINE_WF( fastq_files,params.host)
    } 
    else if(params.pipeline == "GSV_1") {
        GSV_1_WF( fastq_files)
    }    
    else if(params.pipeline == "GSV_2") {
        Channel
            .fromPath( params.merged_reads, glob:true )
            .ifEmpty { error "No FASTQs match: ${params.merged_reads}" }
        
            /* keep sample-ID (sid), read-type (rtype), and the Path itself */
            .map { f ->
                def m = (f.name =~ /(.+?).(extendedFrags|notCombined)\.fastq\.gz$/)
                if( !m ) error "Unrecognised FLASH name: ${f.name}"
                def sid   = m[0][1]
                def rtype = m[0][2]            // extendedFrags | notCombined
                tuple( sid, tuple(rtype, f) )  // --> (sid , (rtype , file))
            }
        
            /* group by sample ID → (sid , [ (rtype,file) , … ]) */
            .groupTuple()
        
            /* split the list into the two files we need */
            .map { sid, list ->
                def merged_fq   = list.find { it[0] == 'extendedFrags' }?.getAt(1)
                def unmerged_fq = list.find { it[0] == 'notCombined'   }?.getAt(1)
                if( !merged_fq || !unmerged_fq )
                    error "Sample ${sid} is missing merged or unmerged FASTQ"
                tuple( sid, merged_fq, unmerged_fq )
            }
            .set { to_dedup_ch }      // ( sid , merged , unmerged )
        GSV_2_WF( to_dedup_ch)
    }     
    else if(params.pipeline == "GSV_3") {
        Channel
            .fromPath( params.merged_reads , glob:true )
            .ifEmpty { error "No FASTQs match: ${params.merged_reads}" }
            .map { Path f ->
                // capture sample ID and read-type
                def m = (f.name =~ /(.+?)_(merged|unmerged)\.dedup\.fastq\.gz$/)
                if( !m ) error "Bad name for deduped reads: ${f.name}"
                tuple( m[0][1], tuple(m[0][2], f) )      // (sid , (type , file))
            }
            .groupTuple()                                // (sid , [ (type,file) , … ])
            .map { sid, list ->
                def merged_fq   = list.find { it[0] == 'merged'   }?.getAt(1)
                def unmerged_fq = list.find { it[0] == 'unmerged' }?.getAt(1)
                if( !merged_fq || !unmerged_fq )
                    error "Sample '${sid}' missing merged or unmerged FASTQ"
                tuple( sid, merged_fq, unmerged_fq )      // final 3-element tuple
            }
            .set { to_host_rm_ch }
        GSV_3_WF( to_host_rm_ch,params.host)
    }      
    else if(params.pipeline == "GSV_4") {
        Channel
          .fromFilePairs( params.merged_reads, glob: true )
          .ifEmpty { error "No FASTQ files match: ${params.merged_reads}" }
          .map { sample_id, files ->
            //
            // files will be e.g.
            //   [ Path(…/S1_test_merged.dedup.fastq.gz),
            //     Path(…/S1_test_unmerged.dedup.fastq.gz) ]
            //
            def merged   = files.find { it.name.contains('merged')   }
            def unmerged = files.find { it.name.contains('unmerged') }
            assert merged && unmerged : "Sample $sample_id missing one of merged/unmerged"
            tuple( sample_id, merged, unmerged )
          }
          .set { kraken_input_ch }
        GSV_4_WF( kraken_input_ch)
    }
    else if(params.pipeline == "GSV_5") {
                
        /*  Build ( sid , mergedFile , unmergedFile )  -------------------------- */
        Channel
          .fromFilePairs( params.merged_reads, glob: true )
          .ifEmpty { error "No FASTQ files match: ${params.merged_reads}" }
          .map { sample_id, files ->
            //
            // files will be e.g.
            //   [ Path(…/S1_test_merged.dedup.fastq.gz),
            //     Path(…/S1_test_unmerged.dedup.fastq.gz) ]
            //
            def merged   = files.find { it.name.contains('merged')   }
            def unmerged = files.find { it.name.contains('unmerged') }
            assert merged && unmerged : "Sample $sample_id missing one of merged/unmerged"
            tuple( sample_id, merged, unmerged )
          }
          .set { themisto_input_ch }
          
        GSV_5_WF( themisto_input_ch)
    }       
    else if(params.pipeline == "merge") {
        FLASH_MERGE_WF( fastq_files)
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
