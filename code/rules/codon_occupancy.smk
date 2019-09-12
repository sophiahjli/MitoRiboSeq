# Use R code to generate codon occupancy tables

rule codon_occpancy:
    input:
        all_codon_count=config["results_dir"] + "/codon_count/All_codoncount_table.txt",
        mito_info=config["mito_info"],
        mito_aminoacid_codon=config["mito_aminoacid_codon"]
    output:
        mito_data_raw_coverage_bygene=config["results_dir"] + "/QC/mito_data_raw_coverage_bygene.csv",
        mito_data_raw_depth_bygene=config["results_dir"] + "/QC/mito_data_raw_depth_bygene.csv",
        mito_data_raw_coverage_summary_bycondition=config["results_dir"] +  "/QC/mito_data_raw_coverage_summary_bycondition.csv",
        mito_data_raw_depth_summary_bycondition=config["results_dir"] + "/QC/mito_data_raw_depth_summary_bycondition.csv",
        figures=directory(config["results_dir"] + "/figures"),
        mito_occupancy_table=config["results_dir"] + "/tables/mito_occupancy_table.csv"
    conda:
        "../envs/r.yml"
    script:
        "../scripts/codon_occpancy.R"

