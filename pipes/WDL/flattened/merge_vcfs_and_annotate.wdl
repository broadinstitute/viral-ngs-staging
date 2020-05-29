version 1.0




workflow merge_vcfs_and_annotate {
    meta {
        description: "Merge VCFs emitted by GATK UnifiedGenotyper and annotate with snpEff."
    }

    input {
        File reference_fasta
    }

    parameter_meta {
        reference_fasta: {
          description: "Reference genome, all segments/chromosomes in one fasta file. Headers must be Genbank accessions.",
          patterns: ["*.fasta"]
        }
    }

    call interhost__merge_vcfs_gatk as merge_vcfs {
        input:
            ref_fasta = reference_fasta
    }
    call intrahost__annotate_vcf_snpeff as annotate_vcf {
        input:
            ref_fasta = reference_fasta,
            in_vcf    = merge_vcfs.merged_vcf_gz
    }
    output {
        File merged_vcf_gz           = merge_vcfs.merged_vcf_gz
        File merged_vcf_gz_tbi       = merge_vcfs.merged_vcf_gz_tbi
        File merged_annot_vcf_gz     = annotate_vcf.annot_vcf_gz
        File merged_annot_vcf_gz_tbi = annotate_vcf.annot_vcf_gz_tbi
        File merged_annot_txt_gz     = annotate_vcf.annot_txt_gz
    }
}



task interhost__merge_vcfs_gatk {
  input {
    Array[File] in_vcfs_gz
    File        ref_fasta

    Int?     machine_mem_gb
    String   docker="quay.io/broadinstitute/viral-phylo:2.1.0.0"

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




task intrahost__annotate_vcf_snpeff {
  input {
    File           in_vcf
    File           ref_fasta

    Array[String]? snpEffRef
    String?        emailAddress

    Int?           machine_mem_gb
    String         docker="quay.io/broadinstitute/viral-phylo:2.1.0.0"

    String         output_basename = basename(basename(in_vcf, ".gz"), ".vcf")
  }

  parameter_meta {
    in_vcf:             { description: "input VCF to annotate with snpEff", patterns: ["*.vcf","*.vcf.gz"]}
    ref_fasta:          { description: "The sequence containing the accession to use for annotation; only used if snpEffRef is not provided.", patterns: ["*.fasta","*.fa"] }
    snpEffRef:          { description: "list of accessions to build/find snpEff database. If this is not provided, the ID from the reference fasta will be used (it must be a GenBank accession)" }
    emailAddress:       { description: "email address passed to NCBI if we need to download reference sequences" }
  }

  command {
    set -ex -o pipefail

    intrahost.py --version | tee VERSION

    providedSnpRefAccessions="${sep=' ' snpEffRef}"
    if [ -n "$providedSnpRefAccessions" ]; then 
      snpRefAccessions="$providedSnpRefAccessions";
    else
      snpRefAccessions="$(python -c "from Bio import SeqIO; print(' '.join(list(s.id for s in SeqIO.parse('${ref_fasta}', 'fasta'))))")"
    fi
    echo "snpRefAccessions: $snpRefAccessions"

    if (file "${in_vcf}" | grep -q "gzip" ) ; then
      echo "${in_vcf} is already compressed"
    else
      echo "${in_vcf} is not compressed; gzipping..."
      bgzip "${in_vcf}"
    fi
    echo "Creating vcf index"
    tabix -p vcf "${in_vcf}"
        
    interhost.py snpEff \
        "${in_vcf}" \
        $snpRefAccessions \
        "${output_basename}.annot.vcf.gz" \
        ${'--emailAddress=' + emailAddress}

    intrahost.py iSNV_table \
        "${output_basename}.annot.vcf.gz" \
        "${output_basename}.annot.txt.gz"

    tabix -p vcf "${output_basename}.annot.vcf.gz"
  }

  output {
    File        annot_vcf_gz      = "${output_basename}.annot.vcf.gz"
    File        annot_vcf_gz_tbi  = "${output_basename}.annot.vcf.gz.tbi"
    File        annot_txt_gz      = "${output_basename}.annot.txt.gz"
    String      viralngs_version  = read_string("VERSION")
  }
  runtime {
    docker: "${docker}"
    memory: select_first([machine_mem_gb, 4]) + " GB"
    dx_instance_type: "mem1_ssd1_v2_x4"
  }
}


