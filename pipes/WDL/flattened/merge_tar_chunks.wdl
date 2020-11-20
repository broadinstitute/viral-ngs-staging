version 1.0



workflow merge_tar_chunks {
    meta {
        description: "Combine multiple tar files (possibly compressed by gzip, bz2, lz4, zstd, etc) into a single tar file. Originally meant for combining streaming upload chunks from a sequencing run."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
    }

    call demux__merge_tarballs as merge_tarballs
    output {
        File combined_tar = merge_tarballs.combined_tar
    }
}



task demux__merge_tarballs {
  input {
    Array[File]+  tar_chunks
    String        out_filename

    Int?          machine_mem_gb
    String        docker="quay.io/broadinstitute/viral-core:2.1.8"
  }

  command {
    set -ex -o pipefail

    if [ -z "$TMPDIR" ]; then
      export TMPDIR=$(pwd)
    fi

    file_utils.py --version | tee VERSION

    file_utils.py merge_tarballs \
      ${out_filename} ${sep=' ' tar_chunks} \
      --loglevel=DEBUG
  }

  output {
    File    combined_tar      = "${out_filename}"
    String  viralngs_version  = read_string("VERSION")
  }

  runtime {
    docker: "${docker}"
    memory: select_first([machine_mem_gb, 7]) + " GB"
    cpu: 16
    disks: "local-disk 2625 LOCAL"
    dx_instance_type: "mem1_ssd2_v2_x16"
    preemptible: 0
  }
}


