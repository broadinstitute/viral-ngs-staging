version 1.0



workflow augur_from_beast_mcc {
    meta {
        description: "Visualize BEAST output with Nextstrain. This workflow converts a BEAST MCC tree (.tree file) into an Auspice v2 json file. See https://nextstrain-augur.readthedocs.io/en/stable/faq/import-beast.html for details."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
    }

    input {
        File    beast_mcc_tree
    }

    parameter_meta {
        beast_mcc_tree: {
          description: "A maximum clade credibility (MCC) tree (.tree file) that is output from a BEAST run.",
          patterns: ["*.tree"]
        }
    }

    call nextstrain__augur_import_beast as augur_import_beast {
        input:
            beast_mcc_tree = beast_mcc_tree
    }
    call nextstrain__export_auspice_json as export_auspice_json {
        input:
            tree            = augur_import_beast.tree_newick,
            node_data_jsons = [augur_import_beast.node_data_json]
    }

    output {
        File  beast_mcc_tree_newick      = augur_import_beast.tree_newick
        File  auspice_input_json         = export_auspice_json.virus_json
    }
}



task nextstrain__augur_import_beast {
    meta {
        description: "Import BEAST tree into files ready for augur export, including a Newick-formatted tree and node data in json format. See both https://nextstrain-augur.readthedocs.io/en/stable/faq/import-beast.html and https://nextstrain-augur.readthedocs.io/en/stable/usage/cli/import.html for details."
    }
    input {
        File    beast_mcc_tree

        Float?  most_recent_tip_date
        String? tip_date_regex
        String? tip_date_format
        String? tip_date_delimiter

        Int?    machine_mem_gb
        String  docker = "nextstrain/base:build-20200608T223413Z"
    }
    String tree_basename = basename(beast_mcc_tree, ".tree")
    command {
        augur version > VERSION
        AUGUR_RECURSION_LIMIT=10000 augur import beast \
            --mcc "~{beast_mcc_tree}" \
            --output-tree "~{tree_basename}.nwk" \
            --output-node-data "~{tree_basename}.json" \
            ~{"--most-recent-tip-date " + most_recent_tip_date} \
            ~{"--tip-date-regex " + tip_date_regex} \
            ~{"--tip-date-format " + tip_date_format} \
            ~{"--tip-date-delimeter " + tip_date_delimiter}
        cat /proc/uptime | cut -f 1 -d ' ' > UPTIME_SEC
        cat /proc/loadavg | cut -f 3 -d ' ' > LOAD_15M
        cat /sys/fs/cgroup/memory/memory.max_usage_in_bytes > MEM_BYTES
    }
    runtime {
        docker: docker
        memory: select_first([machine_mem_gb, 3]) + " GB"
        cpu :   2
        disks:  "local-disk 50 HDD"
        dx_instance_type: "mem1_ssd1_v2_x2"
        preemptible: 1
    }
    output {
        File   tree_newick    = "~{tree_basename}.nwk"
        File   node_data_json = "~{tree_basename}.json"
        Int    max_ram_gb = ceil(read_float("/sys/fs/cgroup/memory/memory.max_usage_in_bytes")/1000000000)
        Int    runtime_sec = ceil(read_float("UPTIME_SEC"))
        Int    cpu_load_15min = ceil(read_float("LOAD_15M"))
        String augur_version = read_string("VERSION")
    }
}




task nextstrain__export_auspice_json {
    meta {
        description: "export augur files to json suitable for auspice visualization. The metadata tsv input is generally required unless the node_data_jsons comprehensively capture all of it. See https://nextstrain-augur.readthedocs.io/en/stable/usage/cli/export.html"
    }
    input {
        File        auspice_config
        File?       sample_metadata
        File        tree
        Array[File] node_data_jsons

        File?          lat_longs_tsv
        File?          colors_tsv
        Array[String]? geo_resolutions
        Array[String]? color_by_metadata
        File?          description_md
        Array[String]? maintainers
        String?        title

        String docker = "nextstrain/base:build-20200608T223413Z"
    }
    String out_basename = basename(basename(tree, ".nwk"), "_refined_tree")
    command {
        augur version > VERSION
        touch exportargs

        # --node-data
        if [ -n "~{sep=' ' node_data_jsons}" ]; then
            echo "--node-data" >> exportargs
            cat "~{write_lines(node_data_jsons)}" >> exportargs
        fi

        # --geo-resolutions
        VALS="~{write_lines(select_first([geo_resolutions, []]))}"
        if [ -n "$(cat $VALS)" ]; then
            echo "--geo-resolutions" >> exportargs;
        fi
        cat $VALS >> exportargs

        # --color-by-metadata
        VALS="~{write_lines(select_first([color_by_metadata, []]))}"
        if [ -n "$(cat $VALS)" ]; then
            echo "--color-by-metadata" >> exportargs;
        fi
        cat $VALS >> exportargs

        # --title
        if [ -n "~{title}" ]; then
            echo "--title" >> exportargs
            echo "~{title}" >> exportargs
        fi

        # --maintainers
        VALS="~{write_lines(select_first([maintainers, []]))}"
        if [ -n "$(cat $VALS)" ]; then
            echo "--maintainers" >> exportargs;
        fi
        cat $VALS >> exportargs

        (export AUGUR_RECURSION_LIMIT=10000; cat exportargs | tr '\n' '\0' | xargs -0 -t augur export v2 \
            --tree ~{tree} \
            ~{"--metadata " + sample_metadata} \
            --auspice-config ~{auspice_config} \
            ~{"--lat-longs " + lat_longs_tsv} \
            ~{"--colors " + colors_tsv} \
            ~{"--description " + description_md} \
            --output ~{out_basename}_auspice.json)
        cat /proc/uptime | cut -f 1 -d ' ' > UPTIME_SEC
        cat /proc/loadavg | cut -f 3 -d ' ' > LOAD_15M
        cat /sys/fs/cgroup/memory/memory.max_usage_in_bytes > MEM_BYTES
    }
    runtime {
        docker: docker
        memory: "3 GB"
        cpu :   2
        disks:  "local-disk 100 HDD"
        dx_instance_type: "mem1_ssd1_v2_x2"
        preemptible: 0
    }
    output {
        File   virus_json = "~{out_basename}_auspice.json"
        Int    max_ram_gb = ceil(read_float("MEM_BYTES")/1000000000)
        Int    runtime_sec = ceil(read_float("UPTIME_SEC"))
        Int    cpu_load_15min = ceil(read_float("LOAD_15M"))
        String augur_version = read_string("VERSION")
    }
}


