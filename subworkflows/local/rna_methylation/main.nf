include { XPORE_DATAPREP } from '../../../modules/local/xpore_dataprep/main'
include { XPORE_DIFFMOD  } from '../../../modules/local/xpore_diffmod/main'

//
// Detection and quantification of RNA methylation events — differential (between conditions).
// Reuses the shared nanopolish eventalign from SIGNAL_ANALYSIS: dataprep per sample, then a
// single xpore diffmod comparing conditions taken from meta.condition.
// Requires a proper design: >=2 conditions with >=3 replicates each. Direct-RNA only.
//
workflow RNA_METHYLATION {
    take:
    ch_eventalign   // channel: [ meta, path(eventalign.txt) ] from SIGNAL_ANALYSIS (meta.condition set)

    main:
    versions = Channel.empty()

    // Per-sample dataprep.
    XPORE_DATAPREP(ch_eventalign)
    versions = versions.mix(XPORE_DATAPREP.out.versions)

    // Collect every sample into a single differential comparison.
    ch_diffmod_in = XPORE_DATAPREP.out.dataprep
        .map { meta, dir -> tuple(meta.condition, meta.id, dir) }
        .toList()
        .map { rows ->
            def samples = rows.collect { cond, id, _dir -> [ cond, id ] }
            def dirs    = rows.collect { _cond, _id, dir -> dir }
            tuple([ id: 'xpore' ], samples, dirs)
        }

    XPORE_DIFFMOD(ch_diffmod_in)
    versions = versions.mix(XPORE_DIFFMOD.out.versions)

    emit:
    diffmod  = XPORE_DIFFMOD.out.diffmod    // channel: [ meta, path(diffmod.table) ]
    majority = XPORE_DIFFMOD.out.majority   // channel: [ meta, path(majority_direction_kmer_diffmod.table) ]
    mqc      = XPORE_DIFFMOD.out.mqc        // channel: path(*_xpore_mqc.tsv) for MultiQC
    versions
}
