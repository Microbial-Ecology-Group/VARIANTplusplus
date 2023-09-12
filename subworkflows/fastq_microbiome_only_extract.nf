// Load modules
include { runkraken_extract ; krakenresults ; dlkraken ;runConfirmationKraken;extractedKrakenResults } from '../modules/Microbiome/kraken2.nf'

workflow FASTQ_KRAKEN_EXTRACT_WF {
    take:
        read_pairs_ch
        krakendb
        confirmation_db

    main:
        def default_db_path = "$baseDir/data/kraken_db/minikraken_8GB_20200312/"
        def db_path = file(default_db_path).exists() ? default_db_path : params.kraken_db
        if (db_path == null) {
            dlkraken()
            runkraken_extract(read_pairs_ch, dlkraken.out,krakendb_inter)
        } else {
            kraken_db_ch = Channel.value(db_path)
            runkraken_extract(read_pairs_ch, kraken_db_ch,krakendb_inter)
        }
}
