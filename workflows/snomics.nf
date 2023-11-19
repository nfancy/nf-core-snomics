/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowSnomics.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input
]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

ch_genome_fasta = params.fasta ? file(params.fasta) : []
ch_gtf = params.gtf ? file(params.gtf) : []

ch_aligner = params.aligner
ch_cellranger_index = params.cellranger_index ? file(params.cellranger_index) : []

ch_ensembl_mapping = params.ensembl_mapping ? file(params.ensembl_mapping) : []


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'
include { GTF_GENE_FILTER   } from '../modules/local/gtf_gene_filter'
include { CELLRANGER_ALIGN  } from "../subworkflows/local/align_cellranger"
include { CELLRANGER_ARC_ALIGN  } from "../subworkflows/local/align_cellranger_arc"
include { CELLBENDER } from "../subworkflows/local/cellbender.nf"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                      } from '../modules/nf-core/fastqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { CAT_FASTQ                   } from '../modules/nf-core/cat/fastq/main'
include { MACS2_CALLPEAK              } from '../modules/nf-core/macs2/callpeak/main' 


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary

workflow SNOMICS {

    ch_versions = Channel.empty()


    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    ch_fastq = INPUT_CHECK (
        ch_input, 
        ch_aligner
    ).reads

    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    //
    // MODULE: Concatenate .fastq
    //
    ch_fastq = CAT_FASTQ (
        ch_fastq
    ).reads

    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions)

    //
    // MODULE: Run FastQC
    //
    FASTQC (
        ch_fastq
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    if(!ch_cellranger_index) {

    ch_filter_gtf = GTF_GENE_FILTER ( ch_genome_fasta, ch_gtf ).gtf

    } else {
        
        ch_filter_gtf = ch_gtf
    }
    
    // Run cellranger pipeline
    if (params.aligner == "cellranger") {
        CELLRANGER_ALIGN(
            ch_genome_fasta,
            ch_filter_gtf,
            ch_cellranger_index,
            ch_fastq
        )
        ch_versions = ch_versions.mix(CELLRANGER_ALIGN.out.ch_versions)
        ch_cellranger_h5 = CELLRANGER_ALIGN.out.outs
        ch_cellranger_h5
            .map { meta, outs -> 
                [ meta: meta, h5: outs.findAll { it.endsWith("raw_feature_bc_matrix.h5") }[0]]
            }
            .set { ch_cellranger_h5 }
    }

    // Run cellranger-arc pipeline
    if (params.aligner == "cellranger_arc") {
        CELLRANGER_ARC_ALIGN(
            ch_genome_fasta,
            ch_filter_gtf,
            ch_cellranger_index,
            ch_fastq
        )
        ch_versions = ch_versions.mix(CELLRANGER_ARC_ALIGN.out.ch_versions)
        ch_cellranger_h5 = CELLRANGER_ARC_ALIGN.out.outs
        ch_cellranger_h5
            .map { meta, outs -> 
                [ meta: meta, h5: outs.findAll { it.endsWith("raw_feature_bc_matrix.h5") }[0]]
            }
            .set { ch_cellranger_h5 }

        // Run MACS2 on cellranger-arc atac output
        ch_bam = CELLRANGER_ARC_ALIGN.out.outs    
        ch_bam
            .map { meta, outs -> 
                [ meta: meta, bam: outs.findAll { it.endsWith("atac_possorted_bam.bam") }[0], control: [] ]
            }
            .set { ch_bam }

        MACS2_CALLPEAK (
            ch_bam,
            params.macs2_gsize // need to make availible to user
        )
        ch_versions = ch_versions.mix(MACS2_CALLPEAK.out.versions)
    }
    
    //
    // MODULE: Run cellbender
    //
    CELLBENDER(
        ch_cellranger_h5,
        ch_ensembl_mapping
    )
    ch_versions = ch_versions.mix(CELLBENDER.out.ch_versions.first())



    CUSTOM_DUMPSOFTWAREVERSIONS (
            ch_versions.unique().collectFile(name: 'collated_versions.yml')
        )



}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log)
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
