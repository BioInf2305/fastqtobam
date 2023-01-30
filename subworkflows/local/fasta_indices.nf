//
// create indices using bwa index and samtools faidx, depending on user-inputs
//

include { BWA_INDEX } from '../../modules/nf-core/bwa/index/main'
include { SAMTOOLS_FAIDX   } from '../../modules/nf-core/samtools/faidx/main'

workflow FASTA_INDICES{
    take:
        fastaF
    main:
        if ( !params.skip_bwa_idx ){
                def meta = [:]
                meta.id = "assembly"
                meta.single_end = false
                fasta_meta_pre = fastaF.map{ file -> tuple(meta, file)  }
                BWA_INDEX( 
                    fasta_meta_pre
                    )
                bwa_idx_meta = BWA_INDEX.out.index
            }
        else{
                def meta =   [:]
                meta.id = "assembly"
                meta.single_end = false
                idx_dir = Channel.fromPath( params.dir_bwa_idx, type: 'dir')
                bwa_idx_meta =  [meta, idx_dir]
            }
        if ( !params.skip_samtools_faidx ){
                def meta_sam = [:]
                meta_sam.id = "assembly"
                meta_sam.single_end = false
                faidx_meta_pre = fastaF.map{ file -> tuple(meta_sam, file)  }
                SAMTOOLS_FAIDX( 
                    faidx_meta_pre
                    )
                fa_idx = SAMTOOLS_FAIDX.out.fai.map{ meta, fai -> fai}
                fa_idx.view()
            }
        else{
                fa_idx = Channel.fromPath( params.samtools_faidx)
            }
        
    
    emit:
        bwa_idx_meta                                     // channel: [ val(meta), [ bwa_idx_dir ] ]
        fa_idx                              // channel: path(fa_idx)
        bwa_idx_version = BWA_INDEX.out.versions? BWA_INDEX.out.versions.first():''
        samtools_faidx_version = SAMTOOLS_FAIDX.out.versions? SAMTOOLS_FAIDX.out.versions.first():''
}
