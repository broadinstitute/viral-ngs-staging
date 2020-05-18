version 1.0



workflow beast_to_auspice {
    meta {
        description: "Visualize BEAST output with Nextstrain. This workflow converts a BEAST MCC tree (.tree file) into an Auspice v2 json file. See https://nextstrain-augur.readthedocs.io/en/stable/faq/import-beast.html for details."
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
        File  node_data_json             = augur_import_beast.node_data_json
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
        String  docker = "nextstrain/base:build-20200506T095107Z"
    }
    String tree_basename = basename(beast_mcc_tree, ".tree")
    command {
        augur import beast \
            --mcc "~{beast_mcc_tree}" \
            --output-tree "~{tree_basename}.nwk" \
            --output-node-data "~{tree_basename}.json" \
            ~{"--most-recent-tip-date " + most_recent_tip_date} \
            ~{"--tip-date-regex " + tip_date_regex} \
            ~{"--tip-date-format " + tip_date_format} \
            ~{"--tip-date-delimeter " + tip_date_delimiter}
    }
    runtime {
        docker: docker
        memory: select_first([machine_mem_gb, 3]) + " GB"
        cpu :   2
        disks:  "local-disk 50 HDD"
        dx_instance_type: "mem1_ssd1_v2_x2"
        preemptible: 2
    }
    output {
        File tree_newick    = "~{tree_basename}.nwk"
        File node_data_json = "~{tree_basename}.json"
    }
}




task nextstrain__export_auspice_json {
    meta {
        description: "export augur files to json suitable for auspice visualization. The metadata tsv input is generally required unless the node_data_jsons comprehensively capture all of it. See https://nextstrain-augur.readthedocs.io/en/stable/usage/cli/export.html"
    }
    input {
        File        auspice_config
        File?       metadata
        File        tree
        Array[File] node_data_jsons

        File?       lat_longs_tsv
        File?       colors_tsv

        Int?   machine_mem_gb
        String docker = "nextstrain/base:build-20200506T095107Z"
    }
    String out_basename = basename(basename(tree, ".nwk"), "_refined_tree")
    command {
        NODE_DATA_FLAG=""
        if [ -n "~{sep=' ' node_data_jsons}" ]; then
          NODE_DATA_FLAG="--node-data "
        fi
        augur export v2 --tree ~{tree} \
            ~{"--metadata " + metadata} \
            $NODE_DATA_FLAG ~{sep=' ' node_data_jsons}\
            --auspice-config ~{auspice_config} \
            ~{"--lat-longs " + lat_longs_tsv} \
            ~{"--colors " + colors_tsv} \
            --output ~{out_basename}_auspice.json
    }
    runtime {
        docker: docker
        memory: select_first([machine_mem_gb, 3]) + " GB"
        cpu :   2
        disks:  "local-disk 100 HDD"
        dx_instance_type: "mem1_ssd1_v2_x2"
        preemptible: 2
    }
    output {
        File virus_json = "~{out_basename}_auspice.json"
    }
}


