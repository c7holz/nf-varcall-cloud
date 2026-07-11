process MULTIQC {
    input:
    path '*'

    output:
    path "multiqc_report.html", emit: report

    script:
    "touch multiqc_report.html"

    stub:
    "touch multiqc_report.html"
}
