process NEXTCLADE_DATASET {
  tag        "Downloading Nextclade Dataset: ${dataset_name}"
  label      "process_medium"
  container  'nextstrain/nextclade:3.21.2'

  input:
  val(dataset_name)

  output:
  val(dataset_name), emit: dataset_name
  path "dataset", emit: dataset
  path "logs/${task.process}/*.log", emit: log
  path "versions.yml", emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  """
    mkdir -p nextclade dataset logs/${task.process}
    log=logs/${task.process}/${task.process}.${workflow.sessionId}.log

    date > \$log
    nextclade --version >> \$log
    nextclade_version=\$(nextclade --version)

    echo "Getting nextclade dataset for ${dataset_name}" | tee -a \$log
    nextclade dataset list | tee -a \$log

    nextclade dataset get --name ${dataset_name} --output-dir dataset

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      nextclade: \$(nextclade --version | awk '{print \$NF}')
      dataset: ${dataset_name}
      tag: \$(grep "tag" dataset/pathogen.json | grep tag | sed 's/\"//g' | sed 's/,//g' | awk '{print \$NF}')
      container: ${task.container}
    END_VERSIONS
  """
}

process NEXTCLADE {
  tag        "Clade: ${dataset_name}"
  label      "process_medium"
  container  'nextstrain/nextclade:3.21.2'

  input:
  file(fasta)
  val(dataset_name)
  path(dataset)

  output:
  path "nextclade/*.csv", emit: nextclade_file
  path "nextclade/*", emit: results
  tuple file("nextclade/nextclade.aligned.fasta"), file("nextclade/nextclade.nwk"), emit: prealigned, optional: true
  path "logs/${task.process}/*.log", emit: log
  path "versions.yml", emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args   = task.ext.args   ?: "${params.nextclade_options}"
  def files  = fasta instanceof List ? fasta.join(" ") : fasta
  def sanitized_dataset = dataset_name.replaceAll(/[^a-zA-Z0-9]/,'_')
  def prefix = task.ext.prefix ?: "nextclade_${sanitized_dataset}"
  """
    mkdir -p nextclade dataset logs/${task.process}
    log=logs/${task.process}/${task.process}.${workflow.sessionId}.log

    date > \$log
    nextclade --version >> \$log
    nextclade_version=\$(nextclade --version)

    for fasta in ${files}
    do
      cat \$fasta >> ultimate_fasta.fasta
    done

    nextclade run ${args} \
      --input-dataset ${dataset} \
      --output-all=nextclade/ \
      --jobs ${task.cpus} \
      ultimate_fasta.fasta \
      | tee -a \$log

    cp ultimate_fasta.fasta nextclade/${prefix}.fasta

    # Rename output files to include dataset name for clarity
    for file in nextclade/*; do
      if [ "\$(basename \"\$file\")" != "${prefix}.fasta" ]; then
        mv "\$file" "nextclade/${prefix}_\$(basename \"\$file\")"
      fi
    done

    if [ -f "dataset/pathogen.json" ]
    then
      tag=\$(grep "tag" dataset/pathogen.json | grep tag | sed 's/\"//g' | sed 's/,//g' | awk '{print \$NF}')
    else
      tag="NA"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      nextclade: \$(nextclade --version | awk '{print \$NF}')
      dataset: ${dataset_name}
      tag: \$tag
      container: ${task.container}
    END_VERSIONS
  """
}
