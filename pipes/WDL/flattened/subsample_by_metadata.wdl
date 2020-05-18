version 1.0



workflow subsample_by_metadata {
    meta {
        description: "Filter and subsample a sequence set. See https://nextstrain-augur.readthedocs.io/en/stable/usage/cli/filter.html"
    }

    call nextstrain__filter_subsample_sequences as filter_subsample_sequences
    output {
        File filtered_fasta = filter_subsample_sequences.filtered_fasta
    }
}



task nextstrain__filter_subsample_sequences {
    meta {
        description: "Filter and subsample a sequence set. See https://nextstrain-augur.readthedocs.io/en/stable/usage/cli/filter.html"
    }
    input {
        File     sequences_fasta
        File     sample_metadata_tsv

        Int?     sequences_per_group
        String?  group_by
        File?    include
        File?    exclude

        Boolean  non_nucleotide=true

        String?  min_date
        String?  max_date
        Int?     min_length
        File?    priority
        Int?     subsample_seed
        String?  exclude_where
        String?  include_where

        Int?     machine_mem_gb
        String   docker = "nextstrain/base:build-20200506T095107Z"
    }
    parameter_meta {
        sequences_fasta: {
          description: "Set of sequences in fasta format to subsample using augur filter. These must represent a single chromosome/segment of a genome only.",
          patterns: ["*.fasta", "*.fa"]
        }
        sample_metadata_tsv: {
          description: "Metadata in tab-separated text format. See https://nextstrain-augur.readthedocs.io/en/stable/faq/metadata.html for details.",
          patterns: ["*.txt", "*.tsv"]
        }
    }
    String in_basename = basename(sequences_fasta, ".fasta")
    command {
        augur version > VERSION
        augur filter \
            --sequences ~{sequences_fasta} \
            --metadata ~{sample_metadata_tsv} \
            ~{"--min-date " + min_date} \
            ~{"--max-date " + max_date} \
            ~{"--min-length " + min_length} \
            ~{true="--non-nucleotide " false=""  non_nucleotide} \
            ~{"--exclude " + exclude} \
            ~{"--include " + include} \
            ~{"--priority " + priority} \
            ~{"--sequences-per-group " + sequences_per_group} \
            ~{"--group-by " + group_by} \
            ~{"--subsample-seed " + subsample_seed} \
            ~{"--exclude-where " + exclude_where} \
            ~{"--include-where " + include_where} \
            --output "~{in_basename}.filtered.fasta"
    }
    runtime {
        docker: docker
        memory: "4 GB"
        cpu :   2
        disks:  "local-disk 375 LOCAL"
        dx_instance_type: "mem1_ssd1_v2_x2"
        preemptible: 1
    }
    output {
        File filtered_fasta = "~{in_basename}.filtered.fasta"
        String augur_version = read_string("VERSION")
    }
}


