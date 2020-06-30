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

        Float?   min_date
        Float?   max_date
        Int?     min_length
        File?    priority
        Int?     subsample_seed
        String?  exclude_where
        String?  include_where

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
            --output "~{out_fname}" | tee STDOUT
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


