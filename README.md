This repo contains files used in the exercises in the reproducible research course.


|- conda/ # additional files and directories created during tutorial
|    |- code/
|        |- run_qc.sh
|
|- git_jupyter_docker/ # jupyter creates notebooks directory; docker modifies (?) Dockerfile; git: copy content to other local directory
|    |- code/
|    |   |- make_supplementary.Rmd
|    |
|    |- Dockerfile # skeleton?
|    |- environment.yml
|    |- config.yml
|    |- Snakefile
|
|- rmarkdown/ # creates results directory with pdf
|   |- code/
|   |   |- make_supplementary.Rmd # skeleton file
|   |
|   |- intermediate/
|       |- multiqc/
|       |   |- summary_stats.txt
|       |
|       |- counts.tsv
|
|- snakemake/
    |- environment.yml # as from Conda tutorial
    |- Snakefile # skeleton file
