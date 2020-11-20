version 1.0



workflow diff_genome_sets {

    input {
        Array[File]   genome_set_one
        Array[File]   genome_set_two
    }

    scatter(sample in zip(genome_set_one, genome_set_two)) {
        call reports__compare_two_genomes as compare_two_genomes {
            input:
                genome_one = sample.left,
                genome_two = sample.right,
                out_basename = basename(sample.left, '.fasta')
        }
    }

    call reports__tsv_stack as tsv_stack {
        input:
            input_tsvs = compare_two_genomes.comparison_table,
            out_basename = "diff_genome_sets"
    }

    output {
        File diff = tsv_stack.out_tsv
    }

}



task reports__compare_two_genomes {
  input {
    File          genome_one
    File          genome_two
    String        out_basename

    String        docker="quay.io/broadinstitute/viral-assemble:2.1.4.0"
  }

  command {
    set -ex -o pipefail
    assembly.py --version | tee VERSION
    assembly.py alignment_summary "${genome_one}" "${genome_two}" --outfileName "${out_basename}.txt" --printCounts --loglevel=DEBUG
    cat /proc/uptime | cut -f 1 -d ' ' > UPTIME_SEC
    cat /proc/loadavg > CPU_LOAD
    cat /sys/fs/cgroup/memory/memory.max_usage_in_bytes > MEM_BYTES
  }

  output {
    File   comparison_table = "${out_basename}.txt"
    Int    max_ram_gb = ceil(read_float("MEM_BYTES")/1000000000)
    Int    runtime_sec = ceil(read_float("UPTIME_SEC"))
    String cpu_load = read_string("CPU_LOAD")
    String viralngs_version = read_string("VERSION")
  }

  runtime {
    memory: "3 GB"
    cpu: 2
    docker: "${docker}"
    disks: "local-disk 50 HDD"
    dx_instance_type: "mem1_ssd1_v2_x2"
    preemptible: 1
  }
}




task reports__tsv_stack {
  input {
    Array[File]+   input_tsvs
    String         out_basename
    String         docker="quay.io/broadinstitute/viral-core:2.1.8"
  }

  command {
    csvstack -t --filenames \
      ${sep=' ' input_tsvs} \
      | tr , '\t' \
      > ${out_basename}.txt
  }

  output {
    File   out_tsv = "${out_basename}.txt"
  }

  runtime {
    memory: "1 GB"
    cpu: 1
    docker: "${docker}"
    disks: "local-disk 50 HDD"
    dx_instance_type: "mem1_ssd1_v2_x2"
  }

}


