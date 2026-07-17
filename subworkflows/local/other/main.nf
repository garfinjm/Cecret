include { FREYJA                        } from '../../../modules/local/freyja'
include { FREYJA_AGGREGATE              } from '../../../modules/local/freyja'
include { FREYJA_PATHOGEN               } from '../../../modules/local/freyja'
include { FREYJA_UPDATE                 } from '../../../modules/local/freyja'
include { NEXTCLADE                     } from '../../../modules/local/nextclade'
include { NEXTCLADE_DATASET             } from '../../../modules/local/nextclade'
include { UNZIP                         } from '../../../modules/local/local' 
include { VADR                          } from '../../../modules/local/vadr'

workflow OTHER {
  take:
  ch_fastas // channel: fasta
  ch_bam // channel: [meta, bam]
  ch_reference_genome // channel: fasta
  ch_input_dataset // channel: zipped file
  ch_script // channel: workflow scripts

  main:
  log.info """

Running specific species or custom pathogen analysis. This workflow performs
sequence validation, clade assignment, and lineage abundance estimation for 
organisms other than the default SARS-CoV-2.

Relevant params and their values:

- 'params.species' : ${params.species}
    - Designates subworkflows
- 'params.nextclade_datasets': ${params.nextclade_datasets}
    - Comma-separated list of datasets to download and run Nextclade against
    - When multiple datasets are specified, sequences will be classified by each dataset in parallel
    - See Nextclade documentation at 
      https://docs.nextstrain.org/projects/nextclade/en/stable/user/datasets.html to see 
      available datasets.
- 'params.freyja_pathogen': ${params.freyja_pathogen}
    - Designate which pathogen to download in FREYJA_UPDATE.
    - Will be ignored if value is set to 'SARS-CoV-2'.
    - See Freyja's documentation at 
      https://andersen-lab.github.io/Freyja/src/usage/update.html to see available 
      pathogens.
- 'params.vadr_reference': ${params.vadr_reference}
    - Designate with reference to use in the VADR process and corresponding container.
    - Will be ignored if value is equal to 'sarscov2'.
    - See available images at https://hub.docker.com/r/staphb/vadr/tags.
- 'params.download_nextclade_dataset' : ${params.download_nextclade_dataset}
    - Will download the datasets according to 'params.nextclade_datasets' 
      in NEXTCLADE_DATASET.
- 'params.predownloaded_nextclade_dataset' : ${params.predownloaded_nextclade_dataset}
    - Allows the user to use an existing nextclade dataset in a zipped directory.
    - See https://github.com/UPHL-BioNGS/Cecret/wiki/Usage#nextclade-datasets for more 
      information.

┏━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
┃ process            ┃ description                                                       ┃
┣━━━━━━━━━━━━━━━━━━━━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┫
┃ VADR               ┃ Validates viral sequences and annotates expected features/errors. ┃
┃ NEXTCLADE_DATASET  ┃ Downloads the requested Nextclade dataset for the target pathogen.┃
┃ NEXTCLADE          ┃ Performs clade assignment, mutation calling, and QC.              ┃
┃ FREYJA_UPDATE      ┃ Updates the Freyja database for the specified pathogen.           ┃
┃ FREYJA             ┃ Estimates lineage abundances from BAM files (e.g., wastewater).   ┃
┃ FREYJA_AGGREGATE   ┃ Aggregates Freyja abundance outputs across multiple samples.      ┃
┗━━━━━━━━━━━━━━━━━━━━┻━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

"""

  // create some empty channels for optional results
  ch_versions    = channel.empty()
  ch_for_summary = channel.empty()
  ch_for_multiqc = channel.empty()

  // run vadr only if vadr is expected to run
  if (params.vadr_reference && params.vadr_reference != 'sarscov2') {
    VADR(ch_fastas.collect())
    ch_versions    = VADR.out.versions
    ch_for_summary = ch_for_summary.mix(VADR.out.vadr_file)
  }

  // run freyja only if freyja is expected to run
  if (params.freyja_pathogen && params.freyja_pathogen != 'SARS-CoV-2') {

    if (params.freyja_update?.toString()?.toBoolean()) {
      FREYJA_UPDATE(params.freyja_pathogen)
      ch_versions = ch_versions.mix(FREYJA_UPDATE.out.versions)

      FREYJA_PATHOGEN(ch_bam.map{it -> tuple(it[0], it[1])}.combine(ch_reference_genome).combine(FREYJA_UPDATE.out.db))
      ch_versions = ch_versions.mix(FREYJA_PATHOGEN.out.versions.first())
      ch_freyja_out = FREYJA_PATHOGEN.out.demix
    } else {
      FREYJA(ch_bam.map{it -> tuple(it[0], it[1])}.combine(ch_reference_genome))
      ch_versions = ch_versions.mix(FREYJA.out.versions.first())
      ch_freyja_out = FREYJA.out.demix
    }

    if (params.freyja_aggregate?.toString()?.toBoolean()) {
      FREYJA_AGGREGATE(ch_freyja_out.collect(), ch_script)
      ch_versions    = ch_versions.mix(FREYJA_AGGREGATE.out.versions)
      ch_for_multiqc = ch_for_multiqc.mix(FREYJA_AGGREGATE.out.for_multiqc)
      ch_for_summary = ch_for_summary.mix(FREYJA_AGGREGATE.out.aggregated_freyja_file)
    }
  }

  // run nextclade with multiple datasets using flattened channel approach
  if (params.nextclade && params.download_nextclade_dataset?.toString()?.toBoolean()) {
    
    // Parse comma-separated dataset list
    def datasets = params.nextclade_datasets.split(',').collect { it.trim() }
    
    // Create channel from dataset names and download each one
    ch_dataset_names = Channel.fromList(datasets)
    
    NEXTCLADE_DATASET(ch_dataset_names)
    
    // Combine fastas with each dataset
    ch_fasta_and_dataset = ch_fastas.combine(NEXTCLADE_DATASET.out.dataset_name).combine(NEXTCLADE_DATASET.out.dataset)
    
    // Run NEXTCLADE for each dataset in parallel
    NEXTCLADE(ch_fasta_and_dataset.map { fasta, dataset_name, dataset -> fasta }, 
              ch_fasta_and_dataset.map { fasta, dataset_name, dataset -> dataset_name },
              ch_fasta_and_dataset.map { fasta, dataset_name, dataset -> dataset })
    
    ch_versions    = ch_versions.mix(NEXTCLADE_DATASET.out.versions)
    ch_versions    = ch_versions.mix(NEXTCLADE.out.versions)
    ch_for_multiqc = ch_for_multiqc.mix(NEXTCLADE.out.nextclade_file)
    ch_for_summary = ch_for_summary.mix(NEXTCLADE.out.nextclade_file)
    
  } else if (params.nextclade && params.predownloaded_nextclade_dataset) {
    // Fallback to predownloaded dataset (original behavior)
    UNZIP(ch_input_dataset)
    ch_dataset  = UNZIP.out.dataset
    ch_versions = ch_versions.mix(UNZIP.out.versions)
    
    NEXTCLADE(ch_fastas.collect(), Channel.value('predownloaded'), ch_dataset)
    ch_versions    = ch_versions.mix(NEXTCLADE.out.versions)
    ch_for_multiqc = ch_for_multiqc.mix(NEXTCLADE.out.nextclade_file)
    ch_for_summary = ch_for_summary.mix(NEXTCLADE.out.nextclade_file)
  }

  emit:
    for_multiqc = ch_for_multiqc
    for_summary = ch_for_summary
    versions    = ch_versions
}
