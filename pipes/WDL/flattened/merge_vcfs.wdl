version 1.0



workflow merge_vcfs {
    meta {
        description: "Merge VCFs from multiple samples using GATK3."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
    }
    call interhost__merge_vcfs_gatk as merge_vcfs_gatk
    output {
        File merged_vcf_gz       = merge_vcfs_gatk.merged_vcf_gz
        File merged_vcf_gz_tbi   = merge_vcfs_gatk.merged_vcf_gz_tbi
    }
}



task interhost__merge_vcfs_gatk {
  input {
    Array[File] in_vcfs_gz
    File        ref_fasta

    Int?     machine_mem_gb
    String   docker="quay.io/broadinstitute/viral-phylo:2.1.4.0"

    String   output_prefix = "merged"
  }

  parameter_meta {
    in_vcfs_gz: {
      description: "VCF files to merged; should be (b)gzipped.",
      patterns: ["*.vcf.gz"] 
    }
    ref_fasta: {
      description: "fasta file of reference genome relative to which the input VCF sites were called",
      patterns: ["*.fasta",".fa"]
    }
  }

  command {

    # tabix index input vcfs (must be gzipped)
    parallel -I ,, \
      "tabix -p vcf ,," \
      ::: "${sep=' ' in_vcfs_gz}"

    # index reference to create .fai and .dict indices
    samtools faidx "${ref_fasta}"
    picard CreateSequenceDictionary R="${ref_fasta}" O=$(basename $(basename "${ref_fasta}" .fasta) .fa).dict

    # store input vcf file paths in file
    for invcf in $(echo "${sep=' ' in_vcfs_gz}"); do 
      echo "$invcf" > input_vcfs.list
    done

    # merge
    gatk3 -T CombineVariants -R "${ref_fasta}" -V input_vcfs.list -o "${output_prefix}.vcf" -genotypeMergeOptions UNIQUIFY
    
    # bgzip output
    bgzip "${output_prefix}.vcf"

    # tabix index the vcf to create .tbi file
    tabix -p vcf "${output_prefix}.vcf.gz"
  }

  output {
    File merged_vcf_gz     = "${output_prefix}.vcf.gz"
    File merged_vcf_gz_tbi = "${output_prefix}.vcf.gz.tbi"
  }

  runtime {
    docker: "${docker}"
    memory: select_first([machine_mem_gb, 3]) + " GB"
    cpu: 2
    dx_instance_type: "mem1_ssd1_v2_x2"
  }
}


