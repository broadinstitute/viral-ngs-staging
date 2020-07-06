version 1.0







workflow classify_single {
    meta {
         description: "Runs raw reads through taxonomic classification (Kraken2), human read depletion (based on Kraken2), de novo assembly (SPAdes), and FASTQC/multiQC of reads."
         author: "Broad Viral Genomics"
         email:  "viral-ngs@broadinstitute.org"
    }

    input {
        File  reads_bam

        File  ncbi_taxdump_tgz

        File  spikein_db
        File  trim_clip_db

        File  kraken2_db_tgz
        File  krona_taxonomy_db_kraken2_tgz
    }

    parameter_meta {
        reads_bam: {
          description: "Reads to classify. May be unmapped or mapped or both, paired-end or single-end.",
          patterns: ["*.bam"]
        }
        spikein_db: {
          description: "ERCC spike-in sequences",
          patterns: ["*.fasta", "*.fasta.gz", "*.fasta.zst"]
        }
        trim_clip_db: {
          description: "Adapter sequences to remove via trimmomatic prior to SPAdes assembly",
          patterns: ["*.fasta", "*.fasta.gz", "*.fasta.zst"]
        }
        kraken2_db_tgz: {
          description: "Pre-built Kraken database tarball containing three files: hash.k2d, opts.k2d, and taxo.k2d.",
          patterns: ["*.tar.gz", "*.tar.lz4", "*.tar.bz2", "*.tar.zst"]
        }
        krona_taxonomy_db_kraken2_tgz: {
          description: "Krona taxonomy database containing a single file: taxonomy.tab, or possibly just a compressed taxonomy.tab",
          patterns: ["*.tab.zst", "*.tab.gz", "*.tab", "*.tar.gz", "*.tar.lz4", "*.tar.bz2", "*.tar.zst"]
        }
        ncbi_taxdump_tgz: {
          description: "An NCBI taxdump.tar.gz file that contains, at the minimum, a nodes.dmp and names.dmp file.",
          patterns: ["*.tar.gz", "*.tar.lz4", "*.tar.bz2", "*.tar.zst"]
        }
    }

    call reports__fastqc as fastqc_raw {
        input: reads_bam = reads_bam
    }
    call reports__align_and_count as spikein {
        input:
            reads_bam = reads_bam,
            ref_db = spikein_db
    }
    call metagenomics__kraken2 as kraken2 {
        input:
            reads_bam = reads_bam,
            kraken2_db_tgz = kraken2_db_tgz,
            krona_taxonomy_db_tgz = krona_taxonomy_db_kraken2_tgz
    }
    call metagenomics__filter_bam_to_taxa as deplete {
        input:
            classified_bam = reads_bam,
            classified_reads_txt_gz = kraken2.kraken2_reads_report,
            ncbi_taxonomy_db_tgz = ncbi_taxdump_tgz,
            exclude_taxa = true,
            taxonomic_names = ["Vertebrata"],
            out_filename_suffix = "hs_depleted"
    }
    call reports__fastqc as fastqc_cleaned {
        input: reads_bam = deplete.bam_filtered_to_taxa
    }
    call metagenomics__filter_bam_to_taxa as filter_acellular {
        input:
            classified_bam = reads_bam,
            classified_reads_txt_gz = kraken2.kraken2_reads_report,
            ncbi_taxonomy_db_tgz = ncbi_taxdump_tgz,
            exclude_taxa = true,
            taxonomic_names = ["Vertebrata", "other sequences", "Bacteria"],
            out_filename_suffix = "acellular"
    }
    call read_utils__rmdup_ubam as rmdup_ubam {
       input:
            reads_unmapped_bam = filter_acellular.bam_filtered_to_taxa
    }
    call assembly__assemble as spades {
        input:
            assembler = "spades",
            reads_unmapped_bam = rmdup_ubam.dedup_bam,
            trim_clip_db = trim_clip_db,
            always_succeed = true
    }

    output {
        File cleaned_reads_unaligned_bam  = deplete.bam_filtered_to_taxa
        File deduplicated_reads_unaligned = rmdup_ubam.dedup_bam
        File contigs_fasta                = spades.contigs_fasta

        Int  read_counts_raw                 = deplete.classified_taxonomic_filter_read_count_pre
        Int  read_counts_depleted            = deplete.classified_taxonomic_filter_read_count_post
        Int  read_counts_dedup               = rmdup_ubam.dedup_read_count_post
        Int  read_counts_prespades_subsample = spades.subsample_read_count

        File kraken2_summary_report = kraken2.kraken2_summary_report
        File kraken2_krona_plot     = kraken2.krona_report_html

        String kraken2_viral_classify_version = kraken2.viralngs_version
        String deplete_viral_classify_version = deplete.viralngs_version
        String spades_viral_assemble_version  = spades.viralngs_version
    }
}



task reports__fastqc {
  input {
    File     reads_bam

    String   docker="quay.io/broadinstitute/viral-core:2.1.8"
  }

  String   reads_basename=basename(reads_bam, ".bam")

  command {
    set -ex -o pipefail
    reports.py --version | tee VERSION
    reports.py fastqc ${reads_bam} ${reads_basename}_fastqc.html --out_zip ${reads_basename}_fastqc.zip
  }

  output {
    File   fastqc_html      = "${reads_basename}_fastqc.html"
    File   fastqc_zip      = "${reads_basename}_fastqc.zip"
    String viralngs_version = read_string("VERSION")
  }

  runtime {
    memory: "2 GB"
    cpu: 1
    docker: "${docker}"
    disks: "local-disk 375 LOCAL"
    dx_instance_type: "mem1_ssd1_v2_x2"
  }
}




task reports__align_and_count {
  input {
    File    reads_bam
    File    ref_db
    Int?    topNHits = 3

    Int?    machine_mem_gb
    String  docker="quay.io/broadinstitute/viral-core:2.1.8"
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




task metagenomics__kraken2 {
  meta {
    description: "Runs Kraken2 classification"
  }

  input {
    File     reads_bam
    File     kraken2_db_tgz         # {database.kdb,taxonomy}
    File     krona_taxonomy_db_tgz  # taxonomy.tab
    Float?   confidence_threshold
    Int?     min_base_qual

    Int?     machine_mem_gb
    String   docker="quay.io/broadinstitute/viral-classify:2.1.4.0"
  }

  parameter_meta {
    reads_bam: {
      description: "Reads or contigs to classify. May be unmapped or mapped or both, paired-end or single-end.",
      patterns: ["*.bam", "*.fasta"]
    }
    kraken2_db_tgz: {
      description: "Pre-built Kraken database tarball containing three files: hash.k2d, opts.k2d, and taxo.k2d.",
      patterns: ["*.tar.gz", "*.tar.lz4", "*.tar.bz2", "*.tar.zst"]
    }
    krona_taxonomy_db_tgz: {
      description: "Krona taxonomy database containing a single file: taxonomy.tab, or possibly just a compressed taxonomy.tab",
      patterns: ["*.tab.zst", "*.tab.gz", "*.tab", "*.tar.gz", "*.tar.lz4", "*.tar.bz2", "*.tar.zst"]
    }
    confidence_threshold: {
      description: "Kraken2 confidence score threshold (0.0-1.0). See https://ccb.jhu.edu/software/kraken2/index.shtml?t=manual#confidence-scoring"
    }
    min_base_qual: {
      description: "Minimum base quality used in classification"
    }
  }

  String out_basename=basename(basename(reads_bam, '.bam'), '.fasta')

  command {
    set -ex -o pipefail

    if [ -z "$TMPDIR" ]; then
      export TMPDIR=$(pwd)
    fi
    DB_DIR=$(mktemp -d --suffix _db)
    mkdir -p $DB_DIR/kraken2 $DB_DIR/krona

    # decompress DB to $DB_DIR
    read_utils.py extract_tarball \
      ${kraken2_db_tgz} $DB_DIR/kraken2 \
      --loglevel=DEBUG
    du -hs $DB_DIR/kraken2

    # unpack krona taxonomy.tab
    if [[ ${krona_taxonomy_db_tgz} == *.tar.* ]]; then
      read_utils.py extract_tarball \
        ${krona_taxonomy_db_tgz} $DB_DIR/krona \
        --loglevel=DEBUG &  # we don't need this until later
    else
      if [[ "${krona_taxonomy_db_tgz}" == *.zst ]]; then
        cat "${krona_taxonomy_db_tgz}" | zstd -d > $DB_DIR/krona/taxonomy.tab &
      elif [[ "${krona_taxonomy_db_tgz}" == *.gz ]]; then
        cat "${krona_taxonomy_db_tgz}" | pigz -dc > $DB_DIR/krona/taxonomy.tab &
      elif [[ "${krona_taxonomy_db_tgz}" == *.bz2 ]]; then
        cat "${krona_taxonomy_db_tgz}" | bzip -dc > $DB_DIR/krona/taxonomy.tab &
      else
        cp "${krona_taxonomy_db_tgz}" $DB_DIR/krona/taxonomy.tab &
      fi
    fi

    metagenomics.py --version | tee VERSION

    if [[ ${reads_bam} == *.bam ]]; then
        metagenomics.py kraken2 \
          $DB_DIR/kraken2 \
          ${reads_bam} \
          --outReads   "${out_basename}".kraken2.reads.txt \
          --outReports "${out_basename}".kraken2.report.txt \
          ${"--confidence " + confidence_threshold} \
          ${"--min_base_qual " + min_base_qual} \
          --loglevel=DEBUG
    else # fasta input file: call kraken2 directly
        kraken2 \
          --db $DB_DIR/kraken2 \
          ${reads_bam} \
          --output "${out_basename}".kraken2.reads.txt \
          --report "${out_basename}".kraken2.report.txt \
          ${"--confidence " + confidence_threshold} \
          ${"--min_base_qual " + min_base_qual}
    fi

    wait # for krona_taxonomy_db_tgz to download and extract
    pigz "${out_basename}".kraken2.reads.txt &

    metagenomics.py krona \
      "${out_basename}".kraken2.report.txt \
      $DB_DIR/krona \
      "${out_basename}".kraken2.krona.html \
      --sample_name "${out_basename}" \
      --noRank --noHits --inputType kraken2 \
      --loglevel=DEBUG

    wait # pigz reads.txt
  }

  output {
    File    kraken2_reads_report   = "${out_basename}.kraken2.reads.txt.gz"
    File    kraken2_summary_report = "${out_basename}.kraken2.report.txt"
    File    krona_report_html      = "${out_basename}.kraken2.krona.html"
    String  viralngs_version       = read_string("VERSION")
  }

  runtime {
    docker: "${docker}"
    memory: select_first([machine_mem_gb, 52]) + " GB"
    cpu: 8
    disks: "local-disk 750 LOCAL"
    dx_instance_type: "mem3_ssd1_v2_x8"
    preemptible: 2
  }
}




task metagenomics__filter_bam_to_taxa {
  input {
    File           classified_bam
    File           classified_reads_txt_gz
    File           ncbi_taxonomy_db_tgz # nodes.dmp names.dmp
    Array[String]? taxonomic_names
    Array[Int]?    taxonomic_ids
    Int?           minimum_hit_groups
    Boolean        withoutChildren=false
    Boolean        exclude_taxa=false
    String         out_filename_suffix = "filtered"

    Int?           machine_mem_gb
    String         docker="quay.io/broadinstitute/viral-classify:2.1.4.0"
  }

  String out_basename = basename(classified_bam, ".bam") + "." + out_filename_suffix

  command {
    set -ex -o pipefail
    if [ -z "$TMPDIR" ]; then
      export TMPDIR=$(pwd)
    fi

    # find 90% memory
    mem_in_mb=$(/opt/viral-ngs/source/docker/calc_mem.py mb 90)

    # decompress taxonomy DB to CWD
    read_utils.py extract_tarball \
      ${ncbi_taxonomy_db_tgz} . \
      --loglevel=DEBUG
    if [ -d "taxonomy" ]; then mv taxonomy/* .; fi

    touch taxfilterargs
    TAXNAMELIST="${write_lines(select_first([taxonomic_names, []]))}"
    if [ -n "$(cat $TAXNAMELIST)" ]; then
      echo "--taxNames" >> taxfilterargs
    fi
    cat $TAXNAMELIST >> taxfilterargs

    TAXIDLIST="${write_lines(select_first([taxonomic_ids, []]))}"
    if [ -n "$(cat $TAXIDLIST)" ]; then
      echo "--taxIDs" >> taxfilterargs
    fi
    cat $TAXIDLIST >> taxfilterargs

    metagenomics.py --version | tee VERSION

    samtools view -c ${classified_bam} | tee classified_taxonomic_filter_read_count_pre &

    cat taxfilterargs | xargs -d '\n' metagenomics.py filter_bam_to_taxa \
      ${classified_bam} \
      ${classified_reads_txt_gz} \
      "${out_basename}.bam" \
      nodes.dmp \
      names.dmp \
      ${true='--exclude' false='' exclude_taxa} \
      ${true='--without-children' false='' withoutChildren} \
      ${'--minimum_hit_groups=' + minimum_hit_groups} \
      --out_count COUNT \
      --JVMmemory "$mem_in_mb"m \
      --loglevel=DEBUG

    samtools view -c "${out_basename}.bam" | tee classified_taxonomic_filter_read_count_post
    wait
  }

  output {
    File    bam_filtered_to_taxa                        = "${out_basename}.bam"
    Int     classified_taxonomic_filter_read_count_pre  = read_int("classified_taxonomic_filter_read_count_pre")
    Int     reads_matching_taxa                         = read_int("COUNT")
    Int     classified_taxonomic_filter_read_count_post = read_int("classified_taxonomic_filter_read_count_post")
    String  viralngs_version                            = read_string("VERSION")
  }

  runtime {
    docker: "${docker}"
    memory: select_first([machine_mem_gb, 26]) + " GB"
    disks: "local-disk 375 LOCAL"
    cpu: 4
    dx_instance_type: "mem3_ssd1_v2_x4"
  }

}




task read_utils__rmdup_ubam {
  meta {
    description: "Perform read deduplication on unaligned reads."
  }

  input {
    File     reads_unmapped_bam
    String   method="mvicuna"

    Int?     machine_mem_gb
    String?  docker="quay.io/broadinstitute/viral-core:2.1.8"
  }

  parameter_meta {
    reads_unmapped_bam: { description: "unaligned reads in BAM format", patterns: ["*.bam"] }
    method:             { description: "mvicuna or cdhit" }
  }

  String reads_basename = basename(reads_unmapped_bam, ".bam")
  
  command {
    set -ex -o pipefail
    mem_in_mb=$(/opt/viral-ngs/source/docker/calc_mem.py mb 90)
    read_utils.py --version | tee VERSION

    read_utils.py rmdup_"${method}"_bam \
      "${reads_unmapped_bam}" \
      "${reads_basename}".dedup.bam \
      --JVMmemory "$mem_in_mb"m \
      --loglevel=DEBUG

    samtools view -c ${reads_basename}.dedup.bam | tee dedup_read_count_post
    reports.py fastqc ${reads_basename}.dedup.bam ${reads_basename}.dedup_fastqc.html --out_zip ${reads_basename}.dedup_fastqc.zip
  }

  output {
    File   dedup_bam             = "${reads_basename}.dedup.bam"
    File   dedup_fastqc          = "${reads_basename}.dedup_fastqc.html"
    File   dedup_fastqc_zip      = "${reads_basename}.dedup_fastqc.zip"
    Int    dedup_read_count_post = read_int("dedup_read_count_post")
    String viralngs_version      = read_string("VERSION")
  }

  runtime {
    docker: "${docker}"
    memory: select_first([machine_mem_gb, 7]) + " GB"
    cpu:    2
    disks:  "local-disk 375 LOCAL"
    dx_instance_type: "mem2_ssd1_v2_x2"
  }
}




task assembly__assemble {
    input {
      File     reads_unmapped_bam
      File     trim_clip_db

      Int?     trinity_n_reads=250000
      Int?     spades_n_reads=10000000
      Int?     spades_min_contig_len=0

      String?  assembler="trinity"  # trinity, spades, or trinity-spades
      Boolean? always_succeed=false

      # do this in two steps in case the input doesn't actually have "taxfilt" in the name
      String   sample_name = basename(basename(reads_unmapped_bam, ".bam"), ".taxfilt")

      Int?     machine_mem_gb
      String   docker="quay.io/broadinstitute/viral-assemble:2.1.4.0"
    }

    command {
        set -ex -o pipefail

        # find 90% memory
        mem_in_mb=$(/opt/viral-ngs/source/docker/calc_mem.py mb 90)
        mem_in_gb=$(/opt/viral-ngs/source/docker/calc_mem.py gb 90)

        assembly.py --version | tee VERSION

        if [[ "${assembler}" == "trinity" ]]; then
          assembly.py assemble_trinity \
            ${reads_unmapped_bam} \
            ${trim_clip_db} \
            ${sample_name}.assembly1-${assembler}.fasta \
            ${'--n_reads=' + trinity_n_reads} \
            ${true='--alwaysSucceed' false="" always_succeed} \
            --JVMmemory "$mem_in_mb"m \
            --outReads=${sample_name}.subsamp.bam \
            --loglevel=DEBUG

        elif [[ "${assembler}" == "spades" ]]; then
          assembly.py assemble_spades \
            ${reads_unmapped_bam} \
            ${trim_clip_db} \
            ${sample_name}.assembly1-${assembler}.fasta \
            ${'--nReads=' + spades_n_reads} \
            ${true="--alwaysSucceed" false="" always_succeed} \
            ${'--minContigLen=' + spades_min_contig_len} \
            --memLimitGb $mem_in_gb \
            --outReads=${sample_name}.subsamp.bam \
            --loglevel=DEBUG

        elif [[ "${assembler}" == "trinity-spades" ]]; then
          assembly.py assemble_trinity \
            ${reads_unmapped_bam} \
            ${trim_clip_db} \
            ${sample_name}.assembly1-trinity.fasta \
            ${'--n_reads=' + trinity_n_reads} \
            --JVMmemory "$mem_in_mb"m \
            --outReads=${sample_name}.subsamp.bam \
            ${true='--always_succeed' false='' always_succeed} \
            --loglevel=DEBUG
          assembly.py assemble_spades \
            ${reads_unmapped_bam} \
            ${trim_clip_db} \
            ${sample_name}.assembly1-${assembler}.fasta \
            --contigsUntrusted=${sample_name}.assembly1-trinity.fasta \
            ${'--nReads=' + spades_n_reads} \
            ${true='--alwaysSucceed' false='' always_succeed} \
            ${'--minContigLen=' + spades_min_contig_len} \
            --memLimitGb $mem_in_gb \
            --loglevel=DEBUG

        else
          echo "unrecognized assembler ${assembler}" >&2
          exit 1
        fi

        samtools view -c ${sample_name}.subsamp.bam | tee subsample_read_count >&2
    }

    output {
        File   contigs_fasta        = "${sample_name}.assembly1-${assembler}.fasta"
        File   subsampBam           = "${sample_name}.subsamp.bam"
        Int    subsample_read_count = read_int("subsample_read_count")
        String viralngs_version     = read_string("VERSION")
    }

    runtime {
        docker: "${docker}"
        memory: select_first([machine_mem_gb, 15]) + " GB"
        cpu: 4
        disks: "local-disk 375 LOCAL"
        dx_instance_type: "mem1_ssd1_v2_x8"
    }

}


