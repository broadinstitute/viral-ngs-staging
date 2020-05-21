version 1.0



workflow fetch_annotations {
    call ncbi__download_annotations as download_annotations
    output {
        File        combined_fasta      = download_annotations.combined_fasta
        Array[File] genomes_fasta       = download_annotations.genomes_fasta
        Array[File] features_tbl        = download_annotations.features_tbl
        String      viral_phylo_version = download_annotations.viralngs_version
    }
}



task ncbi__download_annotations {
  input {
    Array[String]+ accessions
    String         emailAddress
    String         combined_out_prefix

    String         docker="quay.io/broadinstitute/viral-phylo:2.0.21.4"
  }

  command {
    set -ex -o pipefail
    ncbi.py --version | tee VERSION
    ncbi.py fetch_feature_tables \
        ${emailAddress} \
        ./ \
        ${sep=' ' accessions} \
        --loglevel DEBUG
    ncbi.py fetch_fastas \
        ${emailAddress} \
        ./ \
        ${sep=' ' accessions} \
        --combinedFilePrefix "${combined_out_prefix}" \
        --loglevel DEBUG
  }

  output {
    File        combined_fasta   = "${combined_out_prefix}.fasta"
    Array[File] genomes_fasta    = glob("*.fasta")
    Array[File] features_tbl     = glob("*.tbl")
    String      viralngs_version = read_string("VERSION")
  }

  runtime {
    docker: "${docker}"
    memory: "7 GB"
    cpu: 2
    dx_instance_type: "mem2_ssd1_v2_x2"
  }
}


