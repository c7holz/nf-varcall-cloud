process MULTIQC {
    input:
    path '*'

    output:
    path "multiqc_report.html", emit: report

    stub:
    "touch multiqc_report.html"
}
