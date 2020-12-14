version 1.0





workflow assemble_refbased {

    meta {
        description: "Reference-based microbial consensus calling. Aligns NGS reads to a singular reference genome, calls a new consensus sequence, and emits: new assembly, reads aligned to provided reference, reads aligned to new assembly, various figures of merit, plots, and QC metrics. The user may provide unaligned reads spread across multiple input files and this workflow will parallelize alignment per input file before merging results prior to consensus calling."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
    }

    parameter_meta {
        sample_name: {
            description: "Base name of output files. The 'SM' field in BAM read group headers are also rewritten to this value. Avoid spaces and other filename-unfriendly characters.",
            category: "common"
        }
        reads_unmapped_bams: {
            description: "Unaligned reads in BAM format",
            patterns: ["*.bam"]
        }
        reference_fasta: {
            description: "Reference genome to align reads to.",
            patterns: ["*.fasta"]
        }
        aligner: {
            description: "Read aligner software to use. Options: novoalign, bwa, minimap2. Minimap2 can automatically handle Illumina, PacBio, or Oxford Nanopore reads as long as the 'PL' field in the BAM read group header is set properly (novoalign and bwa are Illumina-only)."
        }
        novocraft_license: {
            description: "The default Novoalign short read aligner is a commercially licensed software that is available in a much slower, single-threaded version for free. If you have a paid license file, provide it here to run in multi-threaded mode. If this is omitted, it will run in single-threaded mode.",
            patterns: ["*.lic"]
        }
        skip_mark_dupes: {
            description: "skip Picard MarkDuplicates step after alignment. This is recommended to be set to true for PCR amplicon based data. (Default: false)"
        }
        trim_coords_bed: {
            description: "optional primers to trim in reference coordinate space (0-based BED format)",
            patterns: ["*.bed"],
            category: "common"
        }


        assembly_fasta: { description: "The new assembly / consensus sequence for this sample" }
        align_to_ref_variants_vcf_gz: { description: "All variants in the input reads against the original reference genome. This VCF file is used to create the assembly_fasta" }
        assembly_length: { description: "The length of the sequence described in assembly_fasta, inclusive of any uncovered regions denoted by Ns" }
        assembly_length_unambiguous: { description: "The number of called consensus bases in assembly_fasta (excludes regions of the genome that lack read coverage)" }
    }

    input {
        String          sample_name = basename(reads_unmapped_bams[0], '.bam')
        Array[File]+    reads_unmapped_bams
        File            reference_fasta

        String          aligner="minimap2"
        File?           novocraft_license
        Boolean?        skip_mark_dupes=false
        File?           trim_coords_bed
    }

    Map[String,String] align_to_ref_options = {
                            "novoalign": "-r Random -l 40 -g 40 -x 20 -t 501 -k",
                            "bwa": "-k 12 -B 1",
                            "minimap2": ""
                            }
    Map[String,String] align_to_self_options = {
                            "novoalign": "-r Random -l 40 -g 40 -x 20 -t 100",
                            "bwa": "",
                            "minimap2": ""
                            }

    scatter(reads_unmapped_bam in reads_unmapped_bams) {
        call assembly__align_reads as align_to_ref {
            input:
                reference_fasta    = reference_fasta,
                reads_unmapped_bam = reads_unmapped_bam,
                novocraft_license  = novocraft_license,
                skip_mark_dupes    = skip_mark_dupes,
                aligner            = aligner,
                aligner_options    = align_to_ref_options[aligner]
        }
        call assembly__ivar_trim as ivar_trim {
            input:
                aligned_bam = align_to_ref.aligned_only_reads_bam,
                trim_coords_bed = trim_coords_bed
        }
    }

    call read_utils__merge_and_reheader_bams as merge_align_to_ref {
        input:
            in_bams             = ivar_trim.aligned_trimmed_bam,
            sample_name         = sample_name,
            out_basename        = "${sample_name}.align_to_ref.trimmed"
    }

    call assembly__run_discordance as run_discordance {
        input:
            reads_aligned_bam = merge_align_to_ref.out_bam,
            reference_fasta   = reference_fasta,
            out_basename      = sample_name
    }

    call reports__plot_coverage as plot_ref_coverage {
        input:
            aligned_reads_bam   = merge_align_to_ref.out_bam,
            sample_name         = sample_name
    }

    call assembly__refine_assembly_with_aligned_reads as call_consensus {
        input:
            reference_fasta   = reference_fasta,
            reads_aligned_bam = merge_align_to_ref.out_bam,
            sample_name       = sample_name
    }

    scatter(reads_unmapped_bam in reads_unmapped_bams) {
        call assembly__align_reads as align_to_self {
            input:
                reference_fasta    = call_consensus.refined_assembly_fasta,
                reads_unmapped_bam = reads_unmapped_bam,
                novocraft_license  = novocraft_license,
                skip_mark_dupes    = skip_mark_dupes,
                aligner            = aligner,
                aligner_options    = align_to_self_options[aligner]
        }
    }

    call read_utils__merge_and_reheader_bams as merge_align_to_self {
        input:
            in_bams             = align_to_self.aligned_only_reads_bam,
            sample_name         = sample_name,
            out_basename        = "${sample_name}.merge_align_to_self"
    }

    call reports__plot_coverage as plot_self_coverage {
        input:
            aligned_reads_bam   = merge_align_to_self.out_bam,
            sample_name         = sample_name
    }

    output {
        File   assembly_fasta               = call_consensus.refined_assembly_fasta
        File   align_to_ref_variants_vcf_gz = call_consensus.sites_vcf_gz
        Int    assembly_length              = call_consensus.assembly_length
        Int    assembly_length_unambiguous  = call_consensus.assembly_length_unambiguous
        Int    reference_genome_length      = plot_ref_coverage.assembly_length
        Float  assembly_mean_coverage       = plot_ref_coverage.mean_coverage

        Int    dist_to_ref_snps   = call_consensus.dist_to_ref_snps
        Int    dist_to_ref_indels = call_consensus.dist_to_ref_indels

        Int    replicate_concordant_sites  = run_discordance.concordant_sites
        Int    replicate_discordant_snps   = run_discordance.discordant_snps
        Int    replicate_discordant_indels = run_discordance.discordant_indels
        Int    num_read_groups             = run_discordance.num_read_groups
        Int    num_libraries               = run_discordance.num_libraries
        File   replicate_discordant_vcf    = run_discordance.discordant_sites_vcf

        Array[File]   align_to_ref_per_input_aligned_flagstat = align_to_ref.aligned_bam_flagstat
        Array[Int]    align_to_ref_per_input_reads_provided   = align_to_ref.reads_provided
        Array[Int]    align_to_ref_per_input_reads_aligned    = align_to_ref.reads_aligned
        Array[File]   align_to_ref_per_input_fastqc = align_to_ref.aligned_only_reads_fastqc

        File   align_to_ref_merged_aligned_trimmed_only_bam = merge_align_to_ref.out_bam
        File   align_to_ref_merged_coverage_plot            = plot_ref_coverage.coverage_plot
        File   align_to_ref_merged_coverage_tsv             = plot_ref_coverage.coverage_tsv
        Int    align_to_ref_merged_reads_aligned            = plot_ref_coverage.reads_aligned
        Int    align_to_ref_merged_read_pairs_aligned       = plot_ref_coverage.read_pairs_aligned
        Float  align_to_ref_merged_bases_aligned            = plot_ref_coverage.bases_aligned

        File   align_to_self_merged_aligned_only_bam   = merge_align_to_self.out_bam
        File   align_to_self_merged_coverage_plot      = plot_self_coverage.coverage_plot
        File   align_to_self_merged_coverage_tsv       = plot_self_coverage.coverage_tsv
        Int    align_to_self_merged_reads_aligned      = plot_self_coverage.reads_aligned
        Int    align_to_self_merged_read_pairs_aligned = plot_self_coverage.read_pairs_aligned
        Float  align_to_self_merged_bases_aligned      = plot_self_coverage.bases_aligned
        Float  align_to_self_merged_mean_coverage            = plot_self_coverage.mean_coverage

        String align_to_ref_viral_core_version = align_to_ref.viralngs_version[0]
        String ivar_version                    = ivar_trim.ivar_version[0]
        String viral_assemble_version          = call_consensus.viralngs_version
    }

}



task assembly__align_reads {
  meta {
    description: "Align unmapped reads to a reference genome, either using novoalign (default), minimap2, or bwa. Produces an aligned bam file (including all unmapped reads), an aligned-only bam file, both sorted and indexed, along with samtools flagstat output, fastqc stats (on mapped only reads), and some basic figures of merit."
  }

  input {
    File     reference_fasta
    File     reads_unmapped_bam

    File?    novocraft_license

    String   aligner="minimap2"
    String?  aligner_options
    Boolean? skip_mark_dupes=false

    Int?     machine_mem_gb
    String   docker="quay.io/broadinstitute/viral-core:2.1.12"

    String   sample_name = basename(basename(basename(reads_unmapped_bam, ".bam"), ".taxfilt"), ".clean")
  }

  parameter_meta {
    aligner: { description: "Short read aligner to use: novoalign, minimap2, or bwa. (Default: novoalign)" }
  }
  
  command {
    set -ex # do not set pipefail, since grep exits 1 if it can't find the pattern

    read_utils.py --version | tee VERSION

    mem_in_mb=$(/opt/viral-ngs/source/docker/calc_mem.py mb 90)

    cp ${reference_fasta} assembly.fasta
    grep -v '^>' assembly.fasta | tr -d '\n' | wc -c | tee assembly_length

    if [ "$(cat assembly_length)" != "0" ]; then

      # only perform the following if the reference is non-empty

      if [ "${aligner}" == "novoalign" ]; then
        read_utils.py novoindex \
          assembly.fasta \
          ${"--NOVOALIGN_LICENSE_PATH=" + novocraft_license} \
          --loglevel=DEBUG
      fi
      read_utils.py index_fasta_picard assembly.fasta --loglevel=DEBUG
      read_utils.py index_fasta_samtools assembly.fasta --loglevel=DEBUG

      read_utils.py align_and_fix \
        ${reads_unmapped_bam} \
        assembly.fasta \
        --outBamAll "${sample_name}.all.bam" \
        --outBamFiltered "${sample_name}.mapped.bam" \
        --aligner ${aligner} \
        ${'--aligner_options "' + aligner_options + '"'} \
        ${true='--skipMarkDupes' false="" skip_mark_dupes} \
        --JVMmemory "$mem_in_mb"m \
        ${"--NOVOALIGN_LICENSE_PATH=" + novocraft_license} \
        --loglevel=DEBUG

    else
      # handle special case of empty reference fasta -- emit empty bams (with original bam headers)
      samtools view -H -b "${reads_unmapped_bam}" > "${sample_name}.all.bam"
      samtools view -H -b "${reads_unmapped_bam}" > "${sample_name}.mapped.bam"

      samtools index "${sample_name}.all.bam" "${sample_name}.all.bai"
      samtools index "${sample_name}.mapped.bam" "${sample_name}.mapped.bai"
    fi

    cat /proc/loadavg > CPU_LOAD

    # collect figures of merit
    grep -v '^>' assembly.fasta | tr -d '\nNn' | wc -c | tee assembly_length_unambiguous
    samtools view -c ${reads_unmapped_bam} | tee reads_provided
    samtools view -c ${sample_name}.mapped.bam | tee reads_aligned
    # report only primary alignments 260=exclude unaligned reads and secondary mappings
    samtools view -h -F 260 ${sample_name}.all.bam | samtools flagstat - | tee ${sample_name}.all.bam.flagstat.txt
    grep properly ${sample_name}.all.bam.flagstat.txt | cut -f 1 -d ' ' | tee read_pairs_aligned
    samtools view ${sample_name}.mapped.bam | cut -f10 | tr -d '\n' | wc -c | tee bases_aligned
    python -c "print (float("$(cat bases_aligned)")/"$(cat assembly_length_unambiguous)") if "$(cat assembly_length_unambiguous)">0 else print(0)" > mean_coverage

    # fastqc mapped bam
    reports.py fastqc ${sample_name}.mapped.bam ${sample_name}.mapped_fastqc.html --out_zip ${sample_name}.mapped_fastqc.zip

    cat /proc/uptime | cut -f 1 -d ' ' > UPTIME_SEC
    cat /sys/fs/cgroup/memory/memory.max_usage_in_bytes > MEM_BYTES
  }

  output {
    File   aligned_bam                   = "${sample_name}.all.bam"
    File   aligned_bam_idx               = "${sample_name}.all.bai"
    File   aligned_bam_flagstat          = "${sample_name}.all.bam.flagstat.txt"
    File   aligned_only_reads_bam        = "${sample_name}.mapped.bam"
    File   aligned_only_reads_bam_idx    = "${sample_name}.mapped.bai"
    File   aligned_only_reads_fastqc     = "${sample_name}.mapped_fastqc.html"
    File   aligned_only_reads_fastqc_zip = "${sample_name}.mapped_fastqc.zip"
    Int    reads_provided                = read_int("reads_provided")
    Int    reads_aligned                 = read_int("reads_aligned")
    Int    read_pairs_aligned            = read_int("read_pairs_aligned")
    Float  bases_aligned                 = read_float("bases_aligned")
    Float  mean_coverage                 = read_float("mean_coverage")
    Int    max_ram_gb = ceil(read_float("MEM_BYTES")/1000000000)
    Int    runtime_sec = ceil(read_float("UPTIME_SEC"))
    String cpu_load = read_string("CPU_LOAD")
    String viralngs_version              = read_string("VERSION")
  }

  runtime {
    docker: "${docker}"
    memory: select_first([machine_mem_gb, 15]) + " GB"
    cpu: 8
    disks: "local-disk 375 LOCAL"
    dx_instance_type: "mem1_ssd1_v2_x8"
    preemptible: 1
  }
}




task assembly__ivar_trim {
    meta {
      description: "this runs ivar trim on aligned reads, which results in soft-clipping of alignments"
    }

    input {
      File    aligned_bam
      File?   trim_coords_bed
      Int?    min_keep_length
      Int?    sliding_window
      Int?    min_quality=1

      Int?    machine_mem_gb
      String  docker="andersenlabapps/ivar:1.2.2"
    }

    String  bam_basename=basename(aligned_bam, ".bam")

    parameter_meta {
      aligned_bam:     { description: "aligned reads in BAM format", patterns: ["*.bam"] }
      trim_coords_bed: { description: "optional primers to trim in reference coordinate space (0-based BED format)", patterns: ["*.bed"] }
      min_keep_length: { description: "Minimum length of read to retain after trimming (Default: 30)" }
      sliding_window:  { description: "Width of sliding window for quality trimming (Default: 4)" }
      min_quality:     { description: "Minimum quality threshold for sliding window to pass (Default: 20)" }
    }

    command {
        ivar version | head -1 | tee VERSION
        if [ -f "${trim_coords_bed}" ]; then
          ivar trim -e \
            ${'-b ' + trim_coords_bed} \
            ${'-m ' + min_keep_length} \
            ${'-s ' + sliding_window} \
            ${'-q ' + min_quality} \
            -i ${aligned_bam} -p trim
          samtools sort -@ $(nproc) -m 1000M -o ${bam_basename}.trimmed.bam trim.bam
        else
          echo "skipping ivar trim"
          cp "${aligned_bam}" "${bam_basename}.trimmed.bam"
        fi
    }

    output {
        File   aligned_trimmed_bam = "${bam_basename}.trimmed.bam"
        String ivar_version        = read_string("VERSION")
    }

    runtime {
        docker: "${docker}"
        memory: select_first([machine_mem_gb, 7]) + " GB"
        cpu: 4
        disks: "local-disk 375 LOCAL"
        dx_instance_type: "mem1_ssd1_v2_x4"
    }
}




task read_utils__merge_and_reheader_bams {
    meta {
      description: "Merge and/or reheader bam files using a mapping table. This task can modify read group tags in a BAM header field for single BAM files or as part of a BAM merge operation. The output is a single BAM file (given one or more input BAMs) and a three-column tab delimited text table that defines: the field, the old value, and the new value (e.g. LB, old_lib_name, new_lib_name or SM, old_sample_name, new_sample_name)"
    }

    input {
      Array[File]+    in_bams
      String?         sample_name
      File?           reheader_table
      String          out_basename

      String          docker="quay.io/broadinstitute/viral-core:2.1.12"
    }

    command {
        set -ex -o pipefail

        read_utils.py --version | tee VERSION
        mem_in_mb=$(/opt/viral-ngs/source/docker/calc_mem.py mb 90)

        if [ ${length(in_bams)} -gt 1 ]; then
            read_utils.py merge_bams ${sep=' ' in_bams} merged.bam --JVMmemory="$mem_in_mb"m --loglevel DEBUG
        else
            echo "Skipping merge, only one input file"
            cp ${sep=' ' in_bams} merged.bam
        fi    

        # remap all SM values to user specified value
        if [ -n "${sample_name}" ]; then
          # create sample name remapping table based on existing sample names
          samtools view -H merged.bam | perl -n -e'/SM:(\S+)/ && print "SM\t$1\t'"${sample_name}"'\n"' | sort | uniq >> reheader_table.txt
        fi

        # remap arbitrary headers using user specified table
        if [[ -f "${reheader_table}" ]]; then
          cat "${reheader_table}" >> reheader_table.txt
        fi

        # reheader bam file if requested
        if [ -s reheader_table.txt ]; then
          read_utils.py reheader_bam merged.bam reheader_table.txt "${out_basename}.bam" --loglevel DEBUG
        else
          mv merged.bam "${out_basename}.bam"
        fi
    }

    output {
        File   out_bam          = "${out_basename}.bam"
        String viralngs_version = read_string("VERSION")
    }

    runtime {
        docker: "${docker}"
        memory: "3 GB"
        cpu: 2
        disks: "local-disk 750 LOCAL"
        dx_instance_type: "mem1_ssd2_v2_x4"
        preemptible: 0
    }
}




task assembly__run_discordance {
    meta {
      description: "This step evaluates discordance between sequencing runs of the same sample. The input is a merged, aligned BAM file for a single sample. If multiple runs (read groups) exist, we split the aligned reads by read group and separately evaluate consensus calls per read group using bcftools mpileup and call. A VCF is emitted that describes variation between runs."
    }

    input {
      File     reads_aligned_bam
      File     reference_fasta
      String   out_basename = "run"
      Int      min_coverage=4

      String   docker="quay.io/broadinstitute/viral-core:2.1.12"
    }

    command {
        set -ex -o pipefail

        read_utils.py --version | tee VERSION

        # create 2-col table with read group ids in both cols
        python3 <<CODE
        import tools.samtools
        header = tools.samtools.SamtoolsTool().getHeader("${reads_aligned_bam}")
        rgids = [[x[3:] for x in h if x.startswith('ID:')][0] for h in header if h[0]=='@RG']
        n_rgs = len(rgids)
        with open('readgroups.txt', 'wt') as outf:
          for rg in rgids:
            outf.write(rg+'\t'+rg+'\n')
        n_lbs = len(set([[x[3:] for x in h if x.startswith('LB:')][0] for h in header if h[0]=='@RG']))
        with open('num_read_groups', 'wt') as outf:
          outf.write(str(n_rgs)+'\n')
        with open('num_libraries', 'wt') as outf:
          outf.write(str(n_lbs)+'\n')
        CODE

        # bcftools call snps while treating each RG as a separate sample
        bcftools mpileup \
          -G readgroups.txt -d 10000 -a "FORMAT/DP,FORMAT/AD" \
          -q 1 -m 2 -Ou \
          -f "${reference_fasta}" "${reads_aligned_bam}" \
          | bcftools call \
          -P 0 -m --ploidy 1 \
          --threads $(nproc) \
          -Ov -o everything.vcf

        # mask all GT calls when less than 3 reads
        cat everything.vcf | bcftools filter -e "FMT/DP<${min_coverage}" -S . > filtered.vcf
        cat filtered.vcf | bcftools filter -i "MAC>0" > "${out_basename}.discordant.vcf"

        # tally outputs
        bcftools filter -i 'MAC=0' filtered.vcf | bcftools query -f '%POS\n' | wc -l | tee num_concordant
        bcftools filter -i 'TYPE="snp"'  "${out_basename}.discordant.vcf" | bcftools query -f '%POS\n' | wc -l | tee num_discordant_snps
        bcftools filter -i 'TYPE!="snp"' "${out_basename}.discordant.vcf" | bcftools query -f '%POS\n' | wc -l | tee num_discordant_indels
    }

    output {
        File   discordant_sites_vcf = "${out_basename}.discordant.vcf"
        Int    concordant_sites  = read_int("num_concordant")
        Int    discordant_snps   = read_int("num_discordant_snps")
        Int    discordant_indels = read_int("num_discordant_indels")
        Int    num_read_groups   = read_int("num_read_groups")
        Int    num_libraries     = read_int("num_libraries")
        String viralngs_version  = read_string("VERSION")
    }

    runtime {
        docker: "${docker}"
        memory: "3 GB"
        cpu: 2
        disks: "local-disk 100 HDD"
        dx_instance_type: "mem1_ssd1_v2_x2"
        preemptible: 1
    }
}




task reports__plot_coverage {
  input {
    File     aligned_reads_bam
    String   sample_name

    Boolean skip_mark_dupes=false
    Boolean plot_only_non_duplicates=false
    Boolean bin_large_plots=false
    String?  binning_summary_statistic="max" # max or min

    String   docker="quay.io/broadinstitute/viral-core:2.1.12"
  }
  
  command {
    set -ex -o pipefail

    read_utils.py --version | tee VERSION

    samtools view -c ${aligned_reads_bam} | tee reads_aligned
    if [ "$(cat reads_aligned)" != "0" ]; then
      samtools index -@ "$(nproc)" "${aligned_reads_bam}"

      PLOT_DUPE_OPTION=""
      if [[ "${skip_mark_dupes}" != "true" ]]; then
        PLOT_DUPE_OPTION="${true='--plotOnlyNonDuplicates' false="" plot_only_non_duplicates}"
      fi
      
      BINNING_OPTION="${true='--binLargePlots' false="" bin_large_plots}"

      # plot coverage
      reports.py plot_coverage \
        "${aligned_reads_bam}" \
        "${sample_name}.coverage_plot.pdf" \
        --outSummary "${sample_name}.coverage_plot.txt" \
        --plotFormat pdf \
        --plotWidth 1100 \
        --plotHeight 850 \
        --plotDPI 100 \
        $PLOT_DUPE_OPTION \
        $BINNING_OPTION \
        --binningSummaryStatistic ${binning_summary_statistic} \
        --plotTitle "${sample_name} coverage plot" \
        --loglevel=DEBUG

    else
      touch ${sample_name}.coverage_plot.pdf ${sample_name}.coverage_plot.txt
    fi

    # collect figures of merit
    set +o pipefail # grep will exit 1 if it fails to find the pattern
    samtools view -H ${aligned_reads_bam} | perl -n -e'/^@SQ.*LN:(\d+)/ && print "$1\n"' |  python -c "import sys; print(sum(int(x) for x in sys.stdin))" | tee assembly_length
    # report only primary alignments 260=exclude unaligned reads and secondary mappings
    samtools view -h -F 260 ${aligned_reads_bam} | samtools flagstat - | tee ${sample_name}.flagstat.txt
    grep properly ${sample_name}.flagstat.txt | cut -f 1 -d ' ' | tee read_pairs_aligned
    samtools view ${aligned_reads_bam} | cut -f10 | tr -d '\n' | wc -c | tee bases_aligned
    python -c "print (float("$(cat bases_aligned)")/"$(cat assembly_length)") if "$(cat assembly_length)">0 else print(0)" > mean_coverage
  }

  output {
    File   coverage_plot                 = "${sample_name}.coverage_plot.pdf"
    File   coverage_tsv                  = "${sample_name}.coverage_plot.txt"
    Int    assembly_length               = read_int("assembly_length")
    Int    reads_aligned                 = read_int("reads_aligned")
    Int    read_pairs_aligned            = read_int("read_pairs_aligned")
    Float  bases_aligned                 = read_float("bases_aligned")
    Float  mean_coverage                 = read_float("mean_coverage")
    String viralngs_version              = read_string("VERSION")
  }

  runtime {
    docker: "${docker}"
    memory: "7 GB"
    cpu: 2
    disks: "local-disk 375 LOCAL"
    dx_instance_type: "mem1_ssd1_v2_x4"
    preemptible: 1
  }
}




task assembly__refine_assembly_with_aligned_reads {
    meta {
      description: "This step refines/polishes a genome based on short read alignments, producing a new consensus genome. Uses GATK3 Unified Genotyper to produce new consensus. Produces new genome (fasta), variant calls (VCF), and figures of merit."
    }

    input {
      File     reference_fasta
      File     reads_aligned_bam
      String   sample_name

      Boolean? mark_duplicates=false
      Float?   major_cutoff=0.5
      Int?     min_coverage=3

      Int?     machine_mem_gb
      String   docker="quay.io/broadinstitute/viral-assemble:2.1.12.0"
    }

    parameter_meta {
      major_cutoff: {
        description: "If the major allele is present at a frequency higher than this cutoff, we will call an unambiguous base at that position.  If it is equal to or below this cutoff, we will call an ambiguous base representing all possible alleles at that position."
      }
      min_coverage: {
        description: "Minimum read coverage required to call a position unambiguous."
      }
    }

    command {
        set -ex -o pipefail

        # find 90% memory
        mem_in_mb=$(/opt/viral-ngs/source/docker/calc_mem.py mb 90)

        assembly.py --version | tee VERSION

        if [ ${true='true' false='false' mark_duplicates} == "true" ]; then
          read_utils.py mkdup_picard \
            ${reads_aligned_bam} \
            temp_markdup.bam \
            --JVMmemory "$mem_in_mb"m \
            --loglevel=DEBUG
        else
          ln -s ${reads_aligned_bam} temp_markdup.bam
        fi
        samtools index -@ $(nproc) temp_markdup.bam temp_markdup.bai

        ln -s ${reference_fasta} assembly.fasta
        assembly.py refine_assembly \
          assembly.fasta \
          temp_markdup.bam \
          refined.fasta \
          --already_realigned_bam=temp_markdup.bam \
          --outVcf ${sample_name}.sites.vcf.gz \
          --min_coverage ${min_coverage} \
          --major_cutoff ${major_cutoff} \
          --JVMmemory "$mem_in_mb"m \
          --loglevel=DEBUG

        file_utils.py rename_fasta_sequences \
          refined.fasta "${sample_name}.fasta" "${sample_name}"

        # collect variant counts
        bcftools filter -e "FMT/DP<${min_coverage}" -S . "${sample_name}.sites.vcf.gz" -Ou | bcftools filter -i "AC>1" -Ou > "${sample_name}.diffs.vcf"
        bcftools filter -i 'TYPE="snp"'  "${sample_name}.diffs.vcf" | bcftools query -f '%POS\n' | wc -l | tee num_snps
        bcftools filter -i 'TYPE!="snp"' "${sample_name}.diffs.vcf" | bcftools query -f '%POS\n' | wc -l | tee num_indels

        # collect figures of merit
        set +o pipefail # grep will exit 1 if it fails to find the pattern
        grep -v '^>' refined.fasta | tr -d '\n' | wc -c | tee assembly_length
        grep -v '^>' refined.fasta | tr -d '\nNn' | wc -c | tee assembly_length_unambiguous
    }

    output {
        File   refined_assembly_fasta       = "${sample_name}.fasta"
        File   sites_vcf_gz                 = "${sample_name}.sites.vcf.gz"
        Int    assembly_length              = read_int("assembly_length")
        Int    assembly_length_unambiguous  = read_int("assembly_length_unambiguous")
        Int    dist_to_ref_snps             = read_int("num_snps")
        Int    dist_to_ref_indels           = read_int("num_indels")
        String viralngs_version             = read_string("VERSION")
    }

    runtime {
        docker: "${docker}"
        memory: select_first([machine_mem_gb, 7]) + " GB"
        cpu: 8
        disks: "local-disk 375 LOCAL"
        dx_instance_type: "mem1_ssd1_v2_x8"
    }
}


