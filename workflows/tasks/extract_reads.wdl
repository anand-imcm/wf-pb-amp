version 1.0

# extract hifi reads from cluster read info
task extractClusteredHifiReads {
    
    input {
        File hifi_fastq
        File hifi_fastq_index
        File pbaa_read_info
        String file_label
        String docker
    }

    command <<<
        set -euo pipefail

        ln -s ~{hifi_fastq} hifi_reads.fastq

        ln -s ~{hifi_fastq_index} hifi_reads.fastq.fai
        
        cut -d' ' -f1 ~{pbaa_read_info} > ~{file_label}_clustered_holes.txt
        
        seqtk subseq hifi_reads.fastq ~{file_label}_clustered_holes.txt > ~{file_label}_clustered_hifi_reads.fastq
    >>>

    output {
        File clustered_holes = file_label + "_clustered_holes.txt"
        File clustered_hifi_fastq = file_label + "_clustered_hifi_reads.fastq"
    }

    runtime {
        docker: "~{docker}"
        memory: "16G"
        disks: "local-disk 16 HDD"
    }
}