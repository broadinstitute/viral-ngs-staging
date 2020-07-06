version 1.0



workflow beast_gpu {
    meta {
        description: "Runs BEAST (v1) on a GPU instance. Use with care--this can be expensive if run incorrectly."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
    }

    call interhost__beast as beast
    output {
        File        beast_log    = beast.beast_log
        Array[File] trees        = beast.trees
        String      beast_stdout = beast.beast_stdout
    }
}



task interhost__beast {
  input {
    File     beauti_xml

    String   docker="quay.io/broadinstitute/beast-beagle-cuda:1.10.5pre"
  }

  # TO DO: parameterize gpuType and gpuCount

  command {
    set -e
    beast -beagle_info
    nvidia-smi
    bash -c "sleep 60; nvidia-smi" &
    beast \
      -beagle_multipartition off \
      -beagle_GPU -beagle_cuda -beagle_SSE \
      -beagle_double -beagle_scaling always \
      -beagle_order 1,2,3,4 \
      ${beauti_xml}
  }

  output {
    File        beast_log    = glob("*.log")[0]
    Array[File] trees        = glob("*.trees")
    File        beast_stdout = stdout()
  }

  runtime {
    docker: "${docker}"
    memory: "7 GB"
    cpu:    4
    disks: "local-disk 300 HDD"
    bootDiskSizeGb: 50
    gpu:                 true                # dxWDL
    dx_timeout:          "40H"               # dxWDL
    dx_instance_type:    "mem1_ssd1_gpu2_x8" # dxWDL
    acceleratorType:     "nvidia-tesla-k80"  # GCP PAPIv2
    acceleratorCount:    4                   # GCP PAPIv2
    gpuType:             "nvidia-tesla-k80"  # Terra
    gpuCount:            4                   # Terra
    nvidiaDriverVersion: "396.37"
  }
}


