version 1.0




workflow genbank {

    meta {
        description: "Prepare assemblies for Genbank submission. This includes annotation by simple coordinate transfer from Genbank annotations and a multiple alignment. See https://viral-pipelines.readthedocs.io/en/latest/ncbi_submission.html for details."
    }

    input {
        Array[File]+  reference_fastas
        Array[File]+  reference_feature_tables
        Array[File]+  assemblies_fasta

        File          authors_sbt
        File          biosample_attributes
        Int           taxid
        File?         coverage_table
        String?       sequencingTech
        String?       comment
        String?       organism
        String?       molType='cRNA'
    }

    parameter_meta {
        assemblies_fasta: {
          description: "Genomes to prepare for Genbank submission. One file per genome: all segments/chromosomes included in one file. All fasta files must contain exactly the same number of sequences as reference_fasta (which must equal the number of files in reference_annot_tbl).",
          patterns: ["*.fasta"]
        }
        reference_fastas: {
          description: "Reference genome, each segment/chromosome in a separate fasta file, in the exact same count and order as the segments/chromosomes described in genome_fasta. Headers must be Genbank accessions.",
          patterns: ["*.fasta"]
        }
        reference_feature_tables: {
          description: "NCBI Genbank feature table, each segment/chromosome in a separate TBL file, in the exact same count and order as the segments/chromosomes described in genome_fasta and reference_fastas. Accession numbers in the TBL files must correspond exactly to those in reference_fasta.",
          patterns: ["*.tbl"]
        }
        authors_sbt: {
          description: "A genbank submission template file (SBT) with the author list, created at https://submit.ncbi.nlm.nih.gov/genbank/template/submission/",
          patterns: ["*.sbt"]
        }
        biosample_attributes: {
          description: "A post-submission attributes file from NCBI BioSample, which is available at https://submit.ncbi.nlm.nih.gov/subs/ and clicking on 'Download attributes file with BioSample accessions'.",
          patterns: ["*.txt", "*.tsv"]
        }
        taxid: {
          description: "The NCBI taxonomy ID for the species being submitted in this batch (all sequences in this batch must belong to the same taxid). https://www.ncbi.nlm.nih.gov/taxonomy/"
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
          description: "The type of molecule being described. This defaults to 'cRNA' as this pipeline is most commonly used for viral submissions, but any value allowed by the INSDC controlled vocabulary may be used here. Valid values are described at http://www.insdc.org/controlled-vocabulary-moltype-qualifier",
          category: "common"
        }
        comment: {
          description: "Optional comments that can be displayed in the COMMENT section of the Genbank record. This may include any disclaimers about assembly quality or notes about pre-publication availability or requests to discuss pre-publication use with authors."
        }

    }

    call ncbi__biosample_to_genbank as biosample_to_genbank {
        input:
            biosample_attributes = biosample_attributes,
            num_segments = length(reference_fastas),
            taxid = taxid
    }

    scatter(assembly in assemblies_fasta) {
        call ncbi__align_and_annot_transfer_single as annot {
            input:
                genome_fasta = assembly,
                reference_fastas = reference_fastas,
                reference_feature_tables = reference_feature_tables
        }
    }
 
    call ncbi__prepare_genbank as prep_genbank {
        input:
            assemblies_fasta = assemblies_fasta,
            annotations_tbl = flatten(annot.genome_per_chr_tbls),
            authors_sbt = authors_sbt,
            biosampleMap = biosample_to_genbank.biosample_map,
            genbankSourceTable = biosample_to_genbank.genbank_source_modifier_table,
            coverage_table = coverage_table,
            sequencingTech = sequencingTech,
            comment = comment,
            organism = organism,
            molType = molType
    }

    output {
        File submission_zip = prep_genbank.submission_zip
        File archive_zip    = prep_genbank.archive_zip
        File errorSummary   = prep_genbank.errorSummary

        File biosample_map = biosample_to_genbank.biosample_map
        File genbank_source_table = biosample_to_genbank.genbank_source_modifier_table

        Array[File] transferred_annot_tbls = flatten(annot.genome_per_chr_tbls)
        Array[File] genbank_preview_files      = prep_genbank.genbank_preview_files
        Array[File] validation_files           = prep_genbank.validation_files

        String      viral_phylo_version = prep_genbank.viralngs_version
    }

}



task ncbi__biosample_to_genbank {
  meta {
    description: "Prepares two input metadata files for Genbank submission based on a BioSample registration attributes table (attributes.tsv) since all of the necessary values are there. This produces both a Genbank Source Modifier Table and a BioSample ID map file that can be fed into the prepare_genbank task."
  }
  input {
    File  biosample_attributes
    Int   num_segments=1
    Int   taxid

    String  docker="quay.io/broadinstitute/viral-phylo:2.0.21.4"
  }
  String base = basename(biosample_attributes, ".txt")
  command {
    set -ex -o pipefail
    ncbi.py --version | tee VERSION
    ncbi.py biosample_to_genbank \
        "${biosample_attributes}" \
        ${num_segments} \
        ${taxid} \
        "${base}".genbank.src \
        "${base}".biosample.map.txt \
        --loglevel DEBUG
  }
  output {
    File genbank_source_modifier_table = "${base}.genbank.src"
    File biosample_map                 = "${base}.biosample.map.txt"
  }
  runtime {
    docker: "${docker}"
    memory: "1 GB"
    cpu: 1
    dx_instance_type: "mem1_ssd1_v2_x2"
    preemptible: 1
  }
}




task ncbi__align_and_annot_transfer_single {
  meta {
    description: "Given a reference genome annotation in TBL format (e.g. from Genbank or RefSeq) and new genome not in Genbank, produce new annotation files (TBL format with appropriate coordinate conversions) for the new genome. Resulting output can be fed to tbl2asn for Genbank submission."
  }

  input {
    File         genome_fasta
    Array[File]+ reference_fastas
    Array[File]+ reference_feature_tables

    String  docker="quay.io/broadinstitute/viral-phylo:2.0.21.4"
  }

  parameter_meta {
    genome_fasta: {
      description: "New genome, all segments/chromosomes in one fasta file. Must contain the same number of sequences as reference_fasta",
      patterns: ["*.fasta"]
    }
    reference_fastas: {
      description: "Reference genome, each segment/chromosome in a separate fasta file, in the exact same count and order as the segments/chromosomes described in genome_fasta. Headers must be Genbank accessions.",
      patterns: ["*.fasta"]
    }
    reference_feature_tables: {
      description: "NCBI Genbank feature table, each segment/chromosome in a separate TBL file, in the exact same count and order as the segments/chromosomes described in genome_fasta and reference_fastas. Accession numbers in the TBL files must correspond exactly to those in reference_fasta.",
      patterns: ["*.tbl"]
    }
  }

  command {
    set -e
    ncbi.py --version | tee VERSION
    mkdir -p out
    ncbi.py tbl_transfer_multichr \
        "${genome_fasta}" \
        out \
        --ref_fastas ${sep=' ' reference_fastas} \
        --ref_tbls ${sep=' ' reference_feature_tables} \
        --oob_clip \
        --loglevel DEBUG
  }

  output {
    Array[File]+ genome_per_chr_tbls   = glob("out/*.tbl")
    Array[File]+ genome_per_chr_fastas = glob("out/*.fasta")
    String       viralngs_version      = read_string("VERSION")
  }

  runtime {
    docker: "${docker}"
    memory: "15 GB"
    cpu: 4
    dx_instance_type: "mem2_ssd1_v2_x4"
    preemptible: 2
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
    String?      assembly_method
    String?      assembly_method_version

    Int?         machine_mem_gb
    String       docker="quay.io/broadinstitute/viral-phylo:2.0.21.4"
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
    assembly_method: {
      description: "Very short description of the software approach used to assemble the genome. We typically provide a github link here. If this is specified, assembly_method_version should also be specified."
    }
    assembly_method_version: {
      description: "The version of the software used. If this is specified, assembly_method should also be specified."
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
    if [ -n "${assembly_method}" -a -n "${assembly_method_version}" ]; then
      echo "--assembly_method" >> special_args
      echo "${assembly_method}" >> special_args
      echo "--assembly_method_version" >> special_args
      echo "${assembly_method_version}" >> special_args
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


