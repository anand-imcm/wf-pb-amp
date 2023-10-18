version 1.0

# get cluster barcode QC report
task clusterMetrics {
    
    input {
        File clustered_holes
        File lima_report
        File pbaa_read_info
        String file_label
        String docker
    } 

    command <<<
        grep -f <(cut -d/ -f1,2 ~{clustered_holes}) ~{lima_report} > ~{file_label}_clusters.barcode.report
        
        python /home/anand/Documents/aspire-files/data-oxford/terra.bio/wf-pb-amp/scripts/cluster_qc_summary.py \
            --readinfo ~{pbaa_read_info} \
            --demuxreport ~{file_label}_clusters.barcode.report \
            --prefix ~{file_label}

    >>>

    output {
        File clusterQC_report = file_label + "_clusterQC_report.tsv"
    }

    runtime {
        docker: "~{docker}"
    }
}