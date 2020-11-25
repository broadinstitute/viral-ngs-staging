version 1.0



workflow mafft_and_snp {
    meta {
        description: "Align assemblies with mafft and find SNPs with snp-sites."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
    }

    input {
        Array[File]     assembly_fastas
        File            ref_fasta
        Boolean         run_iqtree=false
    }

    parameter_meta {
        assembly_fastas: {
          description: "Set of assembled genomes to align and build trees. These must represent a single chromosome/segment of a genome only. Fastas may be one-sequence-per-individual or a concatenated multi-fasta (unaligned) or a mixture of the two.",
          patterns: ["*.fasta", "*.fa"]
        }
        ref_fasta: {
          description: "A reference assembly (not included in assembly_fastas) to align assembly_fastas against. Typically from NCBI RefSeq or similar.",
          patterns: ["*.fasta", "*.fa"]
        }
    }

    call nextstrain__concatenate as concatenate {
        input:
            infiles     = assembly_fastas,
            output_name = "all_samples_combined_assembly.fasta"
    }
    call nextstrain__mafft_one_chr as mafft {
        input:
            sequences = concatenate.combined,
            ref_fasta = ref_fasta,
            basename  = "all_samples_aligned.fasta"
    }
    call nextstrain__snp_sites as snp_sites {
        input:
            msa_fasta = mafft.aligned_sequences
    }
    if(run_iqtree) {
        call nextstrain__draft_augur_tree as draft_augur_tree {
            input:
                msa_or_vcf = mafft.aligned_sequences
        }
    }

    output {
        File  combined_assemblies = concatenate.combined
        File  multiple_alignment  = mafft.aligned_sequences
        File  unmasked_snps       = snp_sites.snps_vcf
        File? ml_tree             = draft_augur_tree.aligned_tree
    }
}



task nextstrain__concatenate {
    meta {
        description: "This is nothing more than unix cat."
    }
    input {
        Array[File] infiles
        String      output_name
    }
    command {
        cat ~{sep=" " infiles} > "${output_name}"
    }
    runtime {
        docker: "ubuntu"
        memory: "1 GB"
        cpu:    1
        disks: "local-disk 375 LOCAL"
        dx_instance_type: "mem1_ssd1_v2_x2"
    }
    output {
        File combined = "${output_name}"
    }
}




task nextstrain__mafft_one_chr {
    meta {
        description: "Align multiple sequences from FASTA. Only appropriate for closely related (within 99% nucleotide conservation) genomes. See https://mafft.cbrc.jp/alignment/software/closelyrelatedviralgenomes.html"
    }
    input {
        File     sequences
        File?    ref_fasta
        String   basename
        Boolean  remove_reference = false
        Boolean  keep_length = true
        Boolean  large = false
        Boolean  memsavetree = false

        String   docker = "quay.io/broadinstitute/viral-phylo:2.1.10.0"
        Int      mem_size = 60
        Int      cpus = 32
    }
    command {
        set -e -o pipefail
        touch args.txt

        # boolean options
        echo "~{true='--large' false='' large}" >> args.txt
        echo "~{true='--memsavetree' false='' memsavetree}" >> args.txt
        echo "--auto" >> args.txt

        # if ref_fasta is specified, use "closely related" mode
        # see https://mafft.cbrc.jp/alignment/software/closelyrelatedviralgenomes.html
        if [ -f "~{ref_fasta}" ]; then
            echo --addfragments >> args.txt
            echo "~{sequences}" >> args.txt
            echo "~{ref_fasta}" >> args.txt
        else
            echo "~{sequences}" >> args.txt
        fi

        # mafft align to reference in "closely related" mode
        cat args.txt | grep . | xargs -d '\n' mafft --thread -1 \
            ~{true='--keeplength --mapout' false='' keep_length} \
            > msa.fasta

        # remove reference sequence
        python3 <<CODE
        import Bio.SeqIO
        seq_it = Bio.SeqIO.parse('msa.fasta', 'fasta')
        print("dumping " + str(seq_it.__next__().id))
        Bio.SeqIO.write(seq_it, 'msa_drop_one.fasta', 'fasta')
        CODE
        REMOVE_REF="~{true='--remove-reference' false='' remove_reference}"
        if [ -n "$REMOVE_REF" -a -f "~{ref_fasta}" ]; then
            mv msa_drop_one.fasta "~{basename}_aligned.fasta"
        else
            mv msa.fasta "~{basename}_aligned.fasta"
        fi

        # profiling and stats
        cat /proc/uptime | cut -f 1 -d ' ' > UPTIME_SEC
        cat /proc/loadavg > CPU_LOAD
        cat /sys/fs/cgroup/memory/memory.max_usage_in_bytes > MEM_BYTES
    }
    runtime {
        docker: docker
        memory: mem_size + " GB"
        cpu :   cpus
        disks:  "local-disk 100 HDD"
        preemptible: 0
        dx_instance_type: "mem1_ssd1_v2_x36"
    }
    output {
        File   aligned_sequences = "~{basename}_aligned.fasta"
        Int    max_ram_gb = ceil(read_float("MEM_BYTES")/1000000000)
        Int    runtime_sec = ceil(read_float("UPTIME_SEC"))
        String cpu_load = read_string("CPU_LOAD")
    }
}




task nextstrain__snp_sites {
    input {
        File   msa_fasta
        Boolean allow_wildcard_bases=true
        String docker = "quay.io/biocontainers/snp-sites:2.5.1--hed695b0_0"
    }
    String out_basename = basename(msa_fasta, ".fasta")
    command {
        snp-sites -V > VERSION
        snp-sites -v ~{true="" false="-c" allow_wildcard_bases} -o ~{out_basename}.vcf ~{msa_fasta}
    }
    runtime {
        docker: docker
        memory: "1 GB"
        cpu :   1
        disks:  "local-disk 50 HDD"
        preemptible: 0
        dx_instance_type: "mem1_ssd1_v2_x2"
    }
    output {
        File   snps_vcf = "~{out_basename}.vcf"
        String snp_sites_version = read_string("VERSION")
    }
}




task nextstrain__draft_augur_tree {
    meta {
        description: "Build a tree using iqTree. See https://nextstrain-augur.readthedocs.io/en/stable/usage/cli/tree.html"
    }
    input {
        File     msa_or_vcf

        String   method = "iqtree"
        String   substitution_model = "GTR"
        File?    exclude_sites
        File?    vcf_reference
        String?  tree_builder_args

        Int?     cpus
        String   docker = "nextstrain/base:build-20200629T201240Z"
    }
    parameter_meta {
        msa_or_vcf: {
          description: "Set of alignments (fasta format) or variants (vcf format) to construct a tree from using augur tree (iqTree).",
          patterns: ["*.fasta", "*.fa", "*.vcf", "*.vcf.gz"]
        }
    }
    String out_basename = basename(basename(basename(msa_or_vcf, '.gz'), '.vcf'), '.fasta')
    command {
        set -e
        augur version > VERSION
        AUGUR_RECURSION_LIMIT=10000 augur tree --alignment ~{msa_or_vcf} \
            --output ~{out_basename}_~{method}.nwk \
            --method ~{method} \
            --substitution-model ~{default="GTR" substitution_model} \
            ~{"--exclude-sites " + exclude_sites} \
            ~{"--vcf-reference " + vcf_reference} \
            ~{"--tree-builder-args " + tree_builder_args} \
            --nthreads auto
        cat /proc/uptime | cut -f 1 -d ' ' > UPTIME_SEC
        cat /proc/loadavg > CPU_LOAD
        cat /sys/fs/cgroup/memory/memory.max_usage_in_bytes > MEM_BYTES
    }
    runtime {
        docker: docker
        memory: "32 GB"
        cpu:    select_first([cpus, 64])
        disks:  "local-disk 750 LOCAL"
        dx_instance_type: "mem1_ssd1_v2_x36"
        preemptible: 0
    }
    output {
        File   aligned_tree = "~{out_basename}_~{method}.nwk"
        Int    max_ram_gb = ceil(read_float("MEM_BYTES")/1000000000)
        Int    runtime_sec = ceil(read_float("UPTIME_SEC"))
        String cpu_load = read_string("CPU_LOAD")
        String augur_version = read_string("VERSION")
    }
}


