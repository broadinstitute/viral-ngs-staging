version 1.0



workflow augur_from_newick {
    meta {
        description: "Convert a newick formatted phylogenetic tree into a json suitable for auspice visualization. See https://nextstrain-augur.readthedocs.io/en/stable/usage/cli/export.html"
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
    }

    call nextstrain__export_auspice_json as export_auspice_json
    output {
        File auspice_json = export_auspice_json.virus_json
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
    String out_basename = basename(basename(tree, ".nwk"), "_timetree")
    command {
        set -e -o pipefail
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
        cat /proc/loadavg > CPU_LOAD
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
        String cpu_load = read_string("CPU_LOAD")
        String augur_version = read_string("VERSION")
    }
}


