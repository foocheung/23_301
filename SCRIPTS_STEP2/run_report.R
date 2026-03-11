library(quarto)

quarto::quarto_render("v8_comprehensive_scRNAseq_report.qmd", execute_params = list(config_file = "config5.yaml"))

system("mv v8_comprehensive_scRNAseq_report.html analysis_config5.html")


quarto::quarto_render(
  "v8_comprehensive_scRNAseq_report.qmd",
  execute_params = list(config_file = "config5B.yaml")
)

system("mv v8_comprehensive_scRNAseq_report.html analysis_config5B.html")


quarto::quarto_render(
  "v8_comprehensive_scRNAseq_report.qmd",
  execute_params = list(config_file = "config5C.yaml")
)

system("mv v8_comprehensive_scRNAseq_report.html analysis_config5C.html")




#quarto render v2_comprehensive_scRNAseq_report.qmd \
#-P config_file=config2.yaml