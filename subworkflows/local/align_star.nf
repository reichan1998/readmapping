//
// Align short read (HiC and Illumina) data against the genome
//

include { SAMTOOLS_FASTQ    } from '../../modules/nf-core/samtools/fastq/main'
include { BWAMEM2_MEM       } from '../../modules/nf-core/bwamem2/mem/main'
include { SAMTOOLS_MERGE    } from '../../modules/nf-core/samtools/merge/main'
include { SAMTOOLS_SORMADUP } from '../../modules/local/samtools_sormadup'
include { STAR_ALIGN } from '../modules/nf-core/star/align/main'
include { SAMTOOLS_SORT } from '../../modules/nf-core/samtools/sort/main'



workflow ALIGN_STAR {
    take:
    fasta    // channel: [ val(meta), /path/to/fasta ]
    gtf      // channel: [ val(meta), /path/to/gtf ]
    index    // channel: [ val(meta), /path/to/bwamem2/ ]
    reads    // channel: [ val(meta), /path/to/datafile ]
    star_ignore_sjdbgtf // Boolean
    seq_platform // String
    seq_center // String 



    main:
    ch_versions = Channel.empty()

    // Check file types and branch
    reads
    | branch {
        meta, reads ->
            fastq : reads.findAll { it.getName().toLowerCase() =~ /.*f.*\.gz/ }
            cram : true
    }
    | set { ch_reads }


    // Convert from CRAM to FASTQ only if CRAM files were provided as input
    SAMTOOLS_FASTQ ( ch_reads.cram, false )
    ch_versions = ch_versions.mix ( SAMTOOLS_FASTQ.out.versions.first() )
    
    
    SAMTOOLS_FASTQ.out.fastq
    | mix ( ch_reads.fastq )
    | set { ch_reads_fastq }


    // Map reads with STAR
    STAR_ALIGN ( reads, index, gtf, star_ignore_sjdbgtf, seq_platform, seq_center )
    ch_versions       = ch_versions.mix(STAR_ALIGN.out.versions.first())

    SAMTOOLS_SORT ( STAR_ALIGN.out.bam, fasta )
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions.first())

    // Collect all BWAMEM2 output by sample name
    SAMTOOLS_SORT.out.bam
    | map { meta, bam -> [['id': meta.id.split('_')[0..-2].join('_'), 'datatype': meta.datatype], meta.read_count, bam] }
    | groupTuple( by: [0] )
    | map { meta, read_counts, bams -> [meta + [read_count: read_counts.sum()], bams] }
    | branch {
        meta, bams ->
            single_bam: bams.size() == 1
            multi_bams: true
    }
    | set { ch_bams }


    // Merge, but only if there is more than 1 file
    SAMTOOLS_MERGE ( ch_bams.multi_bams, [ [], [] ], [ [], [] ] )
    ch_versions = ch_versions.mix ( SAMTOOLS_MERGE.out.versions.first() )


    SAMTOOLS_MERGE.out.bam
    | mix ( ch_bams.single_bam )
    | set { ch_bam }


    // Mark duplicates
    SAMTOOLS_SORMADUP ( ch_bam, fasta )
    ch_versions = ch_versions.mix ( SAMTOOLS_SORMADUP.out.versions )

    emit:
    bam      = SAMTOOLS_SORMADUP.out.bam     // channel: [ val(meta), /path/to/bam ]
    versions = ch_versions                   // channel: [ versions.yml ]
}
