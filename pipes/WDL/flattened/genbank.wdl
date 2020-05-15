version 1.0




workflow genbank {

    meta {
        description: "Prepare assemblies for Genbank submission. This includes annotation by simple coordinate transfer from Genbank annotations and a multiple alignment. See https://viral-pipelines.readthedocs.io/en/latest/ncbi_submission.html for details."
    }

    input {
        File          reference_fasta
        Array[File]+  reference_annot_tbl
        Array[File]+  assemblies_fasta

        File          authors_sbt
        File?         biosampleMap
        File?         genbankSourceTable
        File?         coverage_table
        String?       sequencingTech
        String?       comment
        String?       organism
        String?       molType
    }

    parameter_meta {
        assemblies_fasta: {
          description: "Genomes to prepare for Genbank submission. One file per genome: all segments/chromosomes included in one file. All fasta files must contain exactly the same number of sequences as reference_fasta (which must equal the number of files in reference_annot_tbl).",
          patterns: ["*.fasta"]
        }
        reference_fasta: {
          description: "Reference genome, all segments/chromosomes in one fasta file. Headers must be Genbank accessions.",
          patterns: ["*.fasta"]
        }
        reference_annot_tbl: {
          description: "NCBI Genbank feature tables, one file for each segment/chromosome described in reference_fasta.",
          patterns: ["*.tbl"]
        }
        authors_sbt: {
          description: "A genbank submission template file (SBT) with the author list, created at https://submit.ncbi.nlm.nih.gov/genbank/template/submission/",
          patterns: ["*.sbt"]
        }
        biosampleMap: {
          description: "A two column tab text file mapping sample IDs (first column) to NCBI BioSample accession numbers (second column). These typically take the format 'SAMN****' and are obtained by registering your samples first at https://submit.ncbi.nlm.nih.gov/",
          patterns: ["*.txt", "*.tsv"],
          category: "common"
        }
        genbankSourceTable: {
          description: "A tab-delimited text file containing requisite metadata for Genbank (a 'source modifier table'). https://www.ncbi.nlm.nih.gov/WebSub/html/help/genbank-source-table.html",
          patterns: ["*.txt", "*.tsv", "*.src"],
          category: "common"
        }
        coverage_table: {
          description: "A two column tab text file mapping sample IDs (first column) to average sequencing coverage (second column, floating point number).",
          patterns: ["*.txt", "*.tsv"],
          category: "common"
        }
        sequencingTech: {
          description: "The type of sequencer used to generate reads. NCBI has a controlled vocabulary for this value which can be found here: https://submit.ncbi.nlm.nih.gov/structcomment/nongenomes/",
          category: "common"
        }
        organism: {
          description: "The scientific name for the organism being submitted. This is typically the species name and should match the name given by the NCBI Taxonomy database. For more info, see: https://www.ncbi.nlm.nih.gov/Sequin/sequin.hlp.html#Organism",
          category: "common"
        }
        molType: {
          description: "The type of molecule being described. This defaults to 'viral cRNA' as this pipeline is most commonly used for viral submissions, but any value allowed by the INSDC controlled vocabulary may be used here. Valid values are described at http://www.insdc.org/controlled-vocabulary-moltype-qualifier",
          category: "common"
        }
        comment: {
          description: "Optional comments that can be displayed in the COMMENT section of the Genbank record. This may include any disclaimers about assembly quality or notes about pre-publication availability or requests to discuss pre-publication use with authors."
        }

    }

    call interhost__multi_align_mafft_ref as mafft {
        input:
            reference_fasta = reference_fasta,
            assemblies_fasta = assemblies_fasta
    }

    scatter(alignment_by_chr in mafft.alignments_by_chr) {
        call ncbi__annot_transfer as annot {
            input:
                multi_aln_fasta = alignment_by_chr,
                reference_fasta = reference_fasta,
                reference_feature_table = reference_annot_tbl
        }
    }
 
    call ncbi__prepare_genbank as prep_genbank {
        input:
            assemblies_fasta = assemblies_fasta,
            annotations_tbl = flatten(annot.transferred_feature_tables),
            authors_sbt = authors_sbt,
            biosampleMap = biosampleMap,
            genbankSourceTable = genbankSourceTable,
            coverage_table = coverage_table,
            sequencingTech = sequencingTech,
            comment = comment,
            organism = organism,
            molType = molType
    }

    output {
        File        submission_zip = prep_genbank.submission_zip
        File        archive_zip    = prep_genbank.archive_zip
        File        errorSummary   = prep_genbank.errorSummary

        Array[File] alignments_by_chr          = mafft.alignments_by_chr
        Array[File] transferred_feature_tables = flatten(annot.transferred_feature_tables)
        Array[File] genbank_preview_files      = prep_genbank.genbank_preview_files
        Array[File] validation_files           = prep_genbank.validation_files

        String      viral_phylo_version = mafft.viralngs_version
    }

}



task interhost__multi_align_mafft_ref {
  input {
    File           reference_fasta
    Array[File]+   assemblies_fasta # fasta files, one per sample, multiple chrs per file okay
    Int?           mafft_maxIters
    Float?         mafft_ep
    Float?         mafft_gapOpeningPenalty

    Int?           machine_mem_gb
    String         docker="quay.io/broadinstitute/viral-phylo:2.0.21.0"
  }

  String           fasta_basename = basename(reference_fasta, '.fasta')

  command {
    interhost.py --version | tee VERSION
    interhost.py multichr_mafft \
      ${reference_fasta} ${sep=' ' assemblies_fasta} \
      . \
      ${'--ep=' + mafft_ep} \
      ${'--gapOpeningPenalty=' + mafft_gapOpeningPenalty} \
      ${'--maxiters=' + mafft_maxIters} \
      --outFilePrefix align_mafft-${fasta_basename} \
      --preservecase \
      --localpair \
      --sampleNameListFile align_mafft-${fasta_basename}-sample_names.txt \
      --loglevel DEBUG
  }

  output {
    #File         sampleNamesFile    = "align_mafft-${fasta_basename}-sample_names.txt"
    Array[File]+  alignments_by_chr  = glob("align_mafft-${fasta_basename}*.fasta")
    String        viralngs_version   = read_string("VERSION")
  }

  runtime {
    docker: "${docker}"
    memory: select_first([machine_mem_gb, 28]) + " GB"
    cpu: 8
    dx_instance_type: "mem2_ssd1_v2_x8"
  }
}




task ncbi__annot_transfer {
  meta {
    description: "Given a reference genome annotation in TBL format (e.g. from Genbank or RefSeq) and a multiple alignment of that reference to other genomes, produce new annotation files (TBL format with appropriate coordinate conversions) for each sequence in the multiple alignment. Resulting output can be fed to tbl2asn for Genbank submission."
  }

  input {
    File         multi_aln_fasta
    File         reference_fasta
    Array[File]+ reference_feature_table

    String  docker="quay.io/broadinstitute/viral-phylo:2.0.21.0"
  }

  parameter_meta {
    multi_aln_fasta: {
      description: "multiple alignment of sample sequences against a reference genome -- for a single chromosome",
      patterns: ["*.fasta"]
    }
    reference_fasta: {
      description: "Reference genome, all segments/chromosomes in one fasta file. Headers must be Genbank accessions.",
      patterns: ["*.fasta"]
    }
    reference_feature_table: {
      description: "NCBI Genbank feature tables, one file for each segment/chromosome described in reference_fasta.",
      patterns: ["*.tbl"]
    }
  }

  command {
    set -e
    ncbi.py --version | tee VERSION
    ncbi.py tbl_transfer_prealigned \
        ${multi_aln_fasta} \
        ${reference_fasta} \
        ${sep=' ' reference_feature_table} \
        . \
        --oob_clip \
        --loglevel DEBUG
  }

  output {
    Array[File] transferred_feature_tables = glob("*.tbl")
    String      viralngs_version           = read_string("VERSION")
  }

  runtime {
    docker: "${docker}"
    memory: "3 GB"
    cpu: 2
    dx_instance_type: "mem1_ssd1_v2_x2"
  }
}




task ncbi__prepare_genbank {
  meta {
    description: "this task runs NCBI's tbl2asn"
  }

  input {
    Array[File]+ assemblies_fasta
    Array[File]  annotations_tbl
    File         authors_sbt
    File?        biosampleMap
    File?        genbankSourceTable
    File?        coverage_table
    String?      sequencingTech
    String?      comment
    String?      organism
    String?      molType

    Int?         machine_mem_gb
    String       docker="quay.io/broadinstitute/viral-phylo:2.0.21.0"
  }

  parameter_meta {
    assemblies_fasta: {
      description: "Assembled genomes. One chromosome/segment per fasta file.",
      patterns: ["*.fasta"]
    }
    annotations_tbl: {
      description: "Gene annotations in TBL format, one per fasta file. Filename basenames must match the assemblies_fasta basenames. These files are typically output from the ncbi.annot_transfer task.",
      patterns: ["*.tbl"]
    }
    authors_sbt: {
      description: "A genbank submission template file (SBT) with the author list, created at https://submit.ncbi.nlm.nih.gov/genbank/template/submission/",
      patterns: ["*.sbt"]
    }
    biosampleMap: {
      description: "A two column tab text file mapping sample IDs (first column) to NCBI BioSample accession numbers (second column). These typically take the format 'SAMN****' and are obtained by registering your samples first at https://submit.ncbi.nlm.nih.gov/",
      patterns: ["*.txt", "*.tsv"]
    }
    genbankSourceTable: {
      description: "A tab-delimited text file containing requisite metadata for Genbank (a 'source modifier table'). https://www.ncbi.nlm.nih.gov/WebSub/html/help/genbank-source-table.html",
      patterns: ["*.txt", "*.tsv"]
    }
    coverage_table: {
      description: "A two column tab text file mapping sample IDs (first column) to average sequencing coverage (second column, floating point number).",
      patterns: ["*.txt", "*.tsv"]
    }
    sequencingTech: {
      description: "The type of sequencer used to generate reads. NCBI has a controlled vocabulary for this value which can be found here: https://submit.ncbi.nlm.nih.gov/structcomment/nongenomes/"
    }
    organism: {
      description: "The scientific name for the organism being submitted. This is typically the species name and should match the name given by the NCBI Taxonomy database. For more info, see: https://www.ncbi.nlm.nih.gov/Sequin/sequin.hlp.html#Organism"
    }
    molType: {
      description: "The type of molecule being described. Any value allowed by the INSDC controlled vocabulary may be used here. Valid values are described at http://www.insdc.org/controlled-vocabulary-moltype-qualifier"
    }
    comment: {
      description: "Optional comments that can be displayed in the COMMENT section of the Genbank record. This may include any disclaimers about assembly quality or notes about pre-publication availability or requests to discuss pre-publication use with authors."
    }

  }

  command {
    set -ex -o pipefail
    ncbi.py --version | tee VERSION
    cp ${sep=' ' annotations_tbl} .

    touch special_args
    if [ -n "${comment}" ]; then
      echo "--comment" >> special_args
      echo "${comment}" >> special_args
    fi
    if [ -n "${sequencingTech}" ]; then
      echo "--sequencing_tech" >> special_args
      echo "${sequencingTech}" >> special_args
    fi
    if [ -n "${organism}" ]; then
      echo "--organism" >> special_args
      echo "${organism}" >> special_args
    fi
    if [ -n "${molType}" ]; then
      echo "--mol_type" >> special_args
      echo "${molType}" >> special_args
    fi
    if [ -n "${coverage_table}" ]; then
      echo -e "sample\taln2self_cov_median" > coverage_table.txt
      cat ${coverage_table} >> coverage_table.txt
      echo "--coverage_table" >> special_args
      echo coverage_table.txt >> special_args
    fi

    cat special_args | xargs -d '\n' ncbi.py prep_genbank_files \
        ${authors_sbt} \
        ${sep=' ' assemblies_fasta} \
        . \
        ${'--biosample_map ' + biosampleMap} \
        ${'--master_source_table ' + genbankSourceTable} \
        --loglevel DEBUG
    zip sequins_only.zip *.sqn
    zip all_files.zip *.sqn *.cmt *.gbf *.src *.fsa *.val
    mv errorsummary.val errorsummary.val.txt # to keep it separate from the glob
  }

  output {
    File        submission_zip           = "sequins_only.zip"
    File        archive_zip              = "all_files.zip"
    Array[File] sequin_files             = glob("*.sqn")
    Array[File] structured_comment_files = glob("*.cmt")
    Array[File] genbank_preview_files    = glob("*.gbf")
    Array[File] source_table_files       = glob("*.src")
    Array[File] fasta_per_chr_files      = glob("*.fsa")
    Array[File] validation_files         = glob("*.val")
    File        errorSummary             = "errorsummary.val.txt"
    String      viralngs_version         = read_string("VERSION")
  }

  runtime {
    docker: "${docker}"
    memory: select_first([machine_mem_gb, 3]) + " GB"
    cpu: 2
    dx_instance_type: "mem1_ssd1_v2_x2"
  }
}


