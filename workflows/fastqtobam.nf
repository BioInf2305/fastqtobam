/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowFastqtobam.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.fasta ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'
include { FASTA_INDICES } from '../subworkflows/local/fasta_indices'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { TRIMGALORE                  } from '../modules/nf-core/trimgalore/main'
include { SICKLE                      } from '../modules/nf-core/sickle/main'
include { FASTQC                      } from '../modules/nf-core/fastqc/main'
include { MULTIQC as MULTIFASTQC      } from '../modules/nf-core/multiqc/main'
include { BWA_MEM                     } from '../modules/nf-core/bwa/mem/main'
include { SAMBAMBA_MARKDUP            } from '../modules/local/sambamba/markdup/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { QUALIMAP_BAMQC              } from '../modules/nf-core/qualimap/bamqc/main'
include { MULTIQC as MULTIBAMQC       } from '../modules/nf-core/multiqc/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow FASTQTOBAM {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    //
    // SUBWORKFLOW: Deals with bwa index and samtools faidx
    //
    
    fastaF = Channel.fromPath(params.fasta)

    FASTA_INDICES(
        fastaF
     )
    
    bwaidx_ver = params.skip_bwa_idx?'':FASTA_INDICES.out.bwa_idx_version

    faidx_ver = params.skip_samtools_faidx?'':FASTA_INDICES.out.samtools_faidx_version

    if ( bwaidx_ver != ''){
        ch_versions = ch_versions.mix(bwaidx_ver)
        }
    if ( faidx_ver != ''){
        ch_version = ch_versions.mix(faidx_ver)
        }

    //
    // MODULE: Run Trimgalore for trimming of adapters
    //
    
    TRIMGALORE(
        INPUT_CHECK.out.reads
    )

    //
    // MODULE: Run sickle for trimming of low-quality reads
    //
    
    SICKLE(
        TRIMGALORE.out.reads.combine([params.qual_type])
    )

    //
    // MODULE: Run fastqc for assessment of read-quality in fastq files
    //

    FASTQC (
        SICKLE.out.paired_trimmed
    )

    ch_versions = ch_versions.mix(FASTQC.out.versions.first())
    ch_versions = ch_versions.mix(TRIMGALORE.out.versions.first())
    ch_versions = ch_versions.mix(SICKLE.out.versions.first())

    //
    // MODULE: BWA_MEN
    //

    SICKLE.out.paired_trimmed
        .combine(FASTA_INDICES.out.bwa_idx_meta)
        .set { mem_params }

    mem_in1 = mem_params.map{ metaR, reads, metaI, idx -> tuple(metaR, reads) }
    mem_in2 = mem_params.map{ metaR, reads, metaI, idx -> tuple(metaI, idx)   }
    

    BWA_MEM(
        mem_in1,
        mem_in2,
        params.sort_bam
    )

    ch_versions = ch_versions.mix(BWA_MEM.out.versions.first())

    /*
    //
    // MODULE: PICARD_MARKDUPLICATES
    //

    BWA_MEM.out.bam
        .combine(fastaF)
        .combine(FASTA_INDICES.out.fa_idx)
        .set{ markduplicates_params }

    markduplicates_in1 = markduplicates_params.map{ metaB, bam, fasF, faIdx -> tuple( metaB, bam )}
    markduplicates_in2 = markduplicates_params.map{ metaB, bam, fasF, faIdx -> fasF }
    markduplicates_in3 = markduplicates_params.map{ metaB, bam, fasF, faIdx -> faIdx }

    PICARD_MARKDUPLICATES(
        markduplicates_in1,
        markduplicates_in2,
        markduplicates_in3
    )
    */

    //
    // MODULE: SAMBAMBA_MARKDUP
    //
    SAMBAMBA_MARKDUP(
        BWA_MEM.out.bam
    )

    ch_versions = ch_versions.mix(SAMBAMBA_MARKDUP.out.versions.first())
    //
    // MODULE: QUALIMAP_BAMQC
    //
    QUALIMAP_BAMQC(
        SAMBAMBA_MARKDUP.out.bam
    )


    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MULTI_FASTQC
    //
    workflow_summary    = WorkflowFastqtobam.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowFastqtobam.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description)
    ch_methods_description = Channel.value(methods_description)

    ch_multifastqc_files = Channel.empty()
    ch_multifastqc_files = ch_multifastqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multifastqc_files = ch_multifastqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multifastqc_files = ch_multifastqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multifastqc_files = ch_multifastqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    ch_multibamqc_files = Channel.empty()
    ch_multibamqc_files = ch_multibamqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multibamqc_files = ch_multibamqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multibamqc_files = ch_multibamqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multibamqc_files = ch_multibamqc_files.mix(QUALIMAP_BAMQC.out.results.collect{it[1]}.ifEmpty([]))

    MULTIFASTQC (
        ch_multifastqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    MULTIBAMQC (
        ch_multibamqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multifastqc_report = MULTIFASTQC.out.report.toList()
    multibamqc_report  = MULTIBAMQC.out.report.toList()

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multifastqc_report, multibamqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
