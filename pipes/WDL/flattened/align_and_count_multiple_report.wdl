version 1.0



workflow align_and_count_multiple_report {
    meta {
        description: "Count the number of times reads map to provided reference sequences. Useful for counting spike-ins, etc."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
    }

    input {
        Array[File]+ reads_unmapped_bams
        File ref_db
    }

    parameter_meta {
        reads_unmapped_bams: {
            description: "Unaligned reads in BAM format",
            patterns: ["*.bam"]
        }
        ref_db: {
            description: "File containing sequences against which reads should me aligned and counted",
            patterns: ["*.fasta","*.fa"]
        }
    }

    scatter(raw_reads in reads_unmapped_bams) {
        call reports__align_and_count as align_and_count {
            input:
                reads_bam = raw_reads,
                ref_db    = ref_db
        }
    }

    call reports__align_and_count_summary as align_and_count_summary {
        input:
            counts_txt = align_and_count.report
    }

    call reports__align_and_count_summary as align_and_count_summary_top_hits {
        input:
            counts_txt    = align_and_count.report_top_hits,
            output_prefix = "count_summary_top_hits"
    }

    output {
        File report               = align_and_count_summary.count_summary
        File report_top_hits      = align_and_count_summary_top_hits.count_summary
        String viral_core_version = align_and_count_summary.viralngs_version
    }
}



task reports__align_and_count {
  input {
    File    reads_bam
    File    ref_db
    Int?    topNHits = 3

    Int?    machine_mem_gb
    String  docker="quay.io/broadinstitute/viral-core:2.1.7"
  }

  String  reads_basename=basename(reads_bam, ".bam")
  String  ref_basename=basename(ref_db, ".fasta")

  command {
    set -ex -o pipefail

    read_utils.py --version | tee VERSION

    ln -s ${reads_bam} ${reads_basename}.bam
    read_utils.py minimap2_idxstats \
      ${reads_basename}.bam \
      ${ref_db} \
      --outStats ${reads_basename}.count.${ref_basename}.txt.unsorted \
      --loglevel=DEBUG

      sort -b -r -n -k3 ${reads_basename}.count.${ref_basename}.txt.unsorted > ${reads_basename}.count.${ref_basename}.txt
      head -n ${topNHits} ${reads_basename}.count.${ref_basename}.txt > ${reads_basename}.count.${ref_basename}.top_${topNHits}_hits.txt
  }

  output {
    File   report           = "${reads_basename}.count.${ref_basename}.txt"
    File   report_top_hits  = "${reads_basename}.count.${ref_basename}.top_${topNHits}_hits.txt"
    String viralngs_version = read_string("VERSION")
  }

  runtime {
    memory: select_first([machine_mem_gb, 15]) + " GB"
    cpu: 4
    docker: "${docker}"
    disks: "local-disk 375 LOCAL"
    dx_instance_type: "mem1_ssd1_v2_x4"
  }
}




task reports__align_and_count_summary {
  input {
    Array[File]+  counts_txt

    String?       output_prefix="count_summary"

    String        docker="quay.io/broadinstitute/viral-core:2.1.7"
  }

  command {
    set -ex -o pipefail

    reports.py --version | tee VERSION
    reports.py aggregate_alignment_counts ${sep=' ' counts_txt} "${output_prefix}".tsv --loglevel=DEBUG
  }

  output {
    File   count_summary    = "${output_prefix}.tsv"
    String viralngs_version = read_string("VERSION")
  }

  runtime {
    memory: "3 GB"
    cpu: 2
    docker: "${docker}"
    disks: "local-disk 50 HDD"
    dx_instance_type: "mem1_ssd1_v2_x2"
  }
}


