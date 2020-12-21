version 1.0



workflow subsample_by_metadata_with_focal {
    meta {
        description: "Filter and subsample a global sequence set with a bias towards a geographic area of interest."
    }

    parameter_meta {
        sample_metadata_tsv: {
            description: "Tab-separated metadata file that contain binning variables and values. Must contain all samples: output will be filtered to the IDs present in this file.",
            patterns: ["*.txt", "*.tsv"]
        }
        sequences_fasta: {
            description: "Sequences in fasta format.",
            patterns: ["*.fasta"]
        }

        focal_variable: {
            description: "The dataset will be bifurcated based on this column header."
        }
        focal_value: {
            description: "The dataset will be bifurcated based whether the focal_variable column matches this value or not. Rows that match this value are considered to be part of the 'focal' set of interest, rows that do not are part of the 'global' set."
        }

        focal_bin_variable: {
            description: "The focal subset of samples will be evenly subsampled across the discrete values of this column header."
        }
        focal_bin_max: {
            description: "The output will contain no more than this number of focal samples from each discrete value in the focal_bin_variable column."
        }

        global_bin_variable: {
            description: "The global subset of samples will be evenly subsampled across the discrete values of this column header."
        }
        global_bin_max: {
            description: "The output will contain no more than this number of global samples from each discrete value in the global_bin_variable column."
        }
    }

    input {
        File    sample_metadata_tsv
        File    sequences_fasta
        File?   priorities

        String  focal_variable = "region"
        String  focal_value = "North America"

        String  focal_bin_variable = "division"
        Int     focal_bin_max = 50

        String  global_bin_variable = "country"
        Int     global_bin_max = 50
    }

    call nextstrain__filter_subsample_sequences as prefilter {
        input:
            sequences_fasta = sequences_fasta,
            sample_metadata_tsv = sample_metadata_tsv
    }

    call nextstrain__filter_subsample_sequences as subsample_focal {
        input:
            sequences_fasta = prefilter.filtered_fasta,
            sample_metadata_tsv = sample_metadata_tsv,
            exclude_where = ["${focal_variable}!=${focal_value}"],
            sequences_per_group = focal_bin_max,
            group_by = focal_bin_variable,
            priority = priorities
    }

    call nextstrain__filter_subsample_sequences as subsample_global {
        input:
            sequences_fasta = prefilter.filtered_fasta,
            sample_metadata_tsv = sample_metadata_tsv,
            exclude_where = ["${focal_variable}=${focal_value}"],
            sequences_per_group = global_bin_max,
            group_by = global_bin_variable,
            priority = priorities
    }

    call nextstrain__concatenate as cat_fasta {
        input:
            infiles = [
                subsample_focal.filtered_fasta, subsample_global.filtered_fasta
            ],
            output_name = "subsampled.fasta"
    }

    call nextstrain__fasta_to_ids as fasta_to_ids {
        input:
            sequences_fasta = cat_fasta.combined
    }

    output {
        File keep_list            = fasta_to_ids.ids_txt
        File subsampled_sequences = cat_fasta.combined
        Int  focal_kept           = subsample_focal.sequences_out
        Int  global_kept          = subsample_global.sequences_out
        Int  sequences_kept       = subsample_focal.sequences_out + subsample_global.sequences_out
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

        Float?   min_date
        Float?   max_date
        Int?     min_length
        File?    priority
        Int?     subsample_seed
        Array[String]?  exclude_where
        Array[String]?  include_where

        String   docker = "nextstrain/base:build-20200629T201240Z"
    }
    parameter_meta {
        sequences_fasta: {
          description: "Set of sequences (unaligned fasta or aligned fasta -- one sequence per genome) or variants (vcf format) to subsample using augur filter.",
          patterns: ["*.fasta", "*.fa", "*.vcf", "*.vcf.gz"]
        }
        sample_metadata_tsv: {
          description: "Metadata in tab-separated text format. See https://nextstrain-augur.readthedocs.io/en/stable/faq/metadata.html for details.",
          patterns: ["*.txt", "*.tsv"]
        }
    }
    String out_fname = sub(sub(basename(sequences_fasta), ".vcf", ".filtered.vcf"), ".fasta$", ".filtered.fasta")
    command {
        set -e
        augur version > VERSION

        touch wherefile
        VALS="~{write_lines(select_first([exclude_where, []]))}"
        if [ -n "$(cat $VALS)" ]; then
            echo "--exclude-where" >> wherefile
            cat $VALS >> wherefile
        fi
        VALS="~{write_lines(select_first([include_where, []]))}"
        if [ -n "$(cat $VALS)" ]; then
            echo "--include-where" >> wherefile
            cat $VALS >> wherefile
        fi

        set -o pipefail
        cat wherefile | tr '\n' '\0' | xargs -0 -t augur filter \
            --sequences "~{sequences_fasta}" \
            --metadata "~{sample_metadata_tsv}" \
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
            --output "~{out_fname}" | tee STDOUT
        set +o pipefail

        #cat ~{sequences_fasta} | grep \> | wc -l > IN_COUNT
        grep "sequences were dropped during filtering" STDOUT | cut -f 1 -d ' ' > DROP_COUNT
        grep "sequences have been written out to" STDOUT | cut -f 1 -d ' ' > OUT_COUNT
        cat /proc/uptime | cut -f 1 -d ' ' > UPTIME_SEC
        cat /proc/loadavg > CPU_LOAD
        cat /sys/fs/cgroup/memory/memory.max_usage_in_bytes > MEM_BYTES
    }
    runtime {
        docker: docker
        memory: "3 GB"
        cpu :   4
        disks:  "local-disk 100 HDD"
        dx_instance_type: "mem1_ssd1_v2_x4"
        preemptible: 1
    }
    output {
        File   filtered_fasta    = out_fname
        String augur_version     = read_string("VERSION")
        Int    sequences_dropped = read_int("DROP_COUNT")
        Int    sequences_out     = read_int("OUT_COUNT")
        Int    max_ram_gb = ceil(read_float("MEM_BYTES")/1000000000)
        Int    runtime_sec = ceil(read_float("UPTIME_SEC"))
        String cpu_load = read_string("CPU_LOAD")
    }
}




task nextstrain__concatenate {
    meta {
        description: "This is nothing more than unix cat."
    }
    input {
        Array[File] infiles
        String      output_name
    }
    command {
        cat ~{sep=" " infiles} > "${output_name}"
    }
    runtime {
        docker: "ubuntu"
        memory: "1 GB"
        cpu:    1
        disks: "local-disk 375 LOCAL"
        dx_instance_type: "mem1_ssd1_v2_x2"
    }
    output {
        File combined = "${output_name}"
    }
}




task nextstrain__fasta_to_ids {
    meta {
        description: "Return the headers only from a fasta file"
    }
    input {
        File sequences_fasta
    }
    String basename = basename(sequences_fasta, ".fasta")
    command {
        cat "~{sequences_fasta}" | grep \> | cut -c 2- > "~{basename}.txt"
    }
    runtime {
        docker: "ubuntu"
        memory: "1 GB"
        cpu:    1
        disks: "local-disk 375 LOCAL"
        dx_instance_type: "mem1_ssd1_v2_x2"
    }
    output {
        File ids_txt = "~{basename}.txt"
    }
}


