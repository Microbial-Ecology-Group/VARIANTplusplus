// Load modules
include { index } from '../modules/Alignment/bwa'
include { bwa_align ; bwa_rm_contaminant_fq ; HostRemovalStats} from '../modules/Alignment/bwa'

//  Host removal for paired reads
workflow FASTQ_RM_HOST_WF {
    take: 
        hostfasta
        read_pairs_ch
    main:
        // Define reference_index variable
        if (params.host_index) {
            reference_index_files = Channel
                .fromPath(params.host_index, glob: true)
                .ifEmpty { error "No files match --host_index '${params.host_index}'" }
                .toList()
                .map { files ->
                    if (files.size() < 6) {
                        error "Expected 6 host index files, found ${files.size()}. Please provide all 6 files, including the host fasta file. Remember to use * in your path."
                    }
                    files.sort()
                }
        } else {
            index(hostfasta)
            reference_index_ch = index.out
        }   
        bwa_rm_contaminant_fq(reference_index_ch, read_pairs_ch )
        HostRemovalStats(bwa_rm_contaminant_fq.out.host_rm_stats.collect())
    emit:
        nonhost_reads = bwa_rm_contaminant_fq.out.nonhost_reads  
}
