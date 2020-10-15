# Use R code to generate codon occupancy tables

rule codon_occupancy:
    input:
        codon_counts=config["results_dir"] + "/codon_count/All_codoncount_table.txt",
        mitochondria_codon_table=config["mito_aminoacid_codon"]
    output:
        expand(config["results_dir"] + "/figures/{sample}_occupancy_plot.pdf", sample=samples.keys()),
        raw_coverage_bygene=config["results_dir"] + "/qc/mito_data_raw_coverage_bygene.csv",
        raw_depth_bygene=config["results_dir"] + "/qc/mito_data_raw_depth_bygene.csv",
        occupancy_table=config["results_dir"] + "/tables/mito_occupancy_table.csv",
        cumsum_table=config["results_dir"] + "/tables/mito_cumsum_table.csv",
        norm_cumsum_plot=config["results_dir"] + "/figures/norm_cumsum_plot.pdf",
    log:
        log_dir + "/codon_occupancy/codon_occupancy.log"
    params:
        figures_directory=config["results_dir"] + "/figures"
    conda:
        "../envs/r.yml"
    script:
        "../scripts/codon_occpancy.R"

rule heatmaps:
    input:
        codon_counts=config["results_dir"] + "/codon_count/All_codoncount_table.txt"
    output:
        heatmap_pdf=config["results_dir"] + "/figures/Heatmap_codoncount_ordercodonfreq.pdf",
        heatmap_scaled_pdf=config["results_dir"] + "/figures/Heatmap_codoncount_ordercodonfreq_scaled.pdf",
        heatmap_scaled_compressed_pdf=config["results_dir"] + "/figures/Heatmap_codoncount_ordercodonfreq_scaled_compressed.pdf",
    log:
        log_dir + "/heatmaps/heatmaps.log"
    conda:
        "../envs/r.yml"
    script:
        "../scripts/heatmap.R"

