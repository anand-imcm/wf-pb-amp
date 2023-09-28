version 1.0

# Create the index of reference genome using pbmm2
task pbmm {
    input {
        File genome_reference
        String docker
    }

    command <<<
        pbmm2 index ~{genome_reference} ref.mmi --preset SUBREAD
    >>>

    output {
        File index = "ref.mmi"
    }

    runtime {
        docker: "~{docker}"
    }

}