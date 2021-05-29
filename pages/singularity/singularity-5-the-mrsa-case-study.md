As a final example we will use `singularity build` to convert the Docker image
of the MRSA project, that we use as a case study in this course, to
a Singularity image:

```bash
singularity build --remote mrsa_proj.sif docker://nbisweden/workshop-reproducible-research
```

This should result in a file called `mrsa_proj.sif`. 

In the Docker image we included the code needed for the workflow in the
`/course` directory of the image. These files are of course also available
in the Singularity image. However, a Singularity image is read-only (unless
using the sandbox feature), and this will be a problem if we try to run the
workflow within the `/course` directory, since the workflow will produce
files and Snakemake will create a `.snakemake` directory. 

Instead, we need to provide the files externally from our host system and
simply use the Singularity image as the environment to execute the workflow
in (i.e. all the software). 

In your current working directory (`singularity/`) the vital MRSA project
files are already available (`Snakefile`, `config.yml`, `code/header.tex` and 
`code/supplementary_material.Rmd`). 

Since Singularity bind mounts the current working directory we can simply
execute the workflow and generate the output files using:

```bash
singularity run --vm-ram 2048 mrsa_proj.sif
```

This executes the default run command, which is 
`snakemake -rp --configfile config.yml` (as defined in the original 
`Dockerfile`). 

> **Note** <br>
> Note here that we have increased the allocated RAM to 2048 MiB (`--vm-ram 2048`), 
> needed to fully run through the workflow. In case the command fails, 
> you can try to increase the RAM to e.g. 4096 MiB, or you can try to run the
> command without the  `--vm-ram` parameter.

The previous step in this tutorial included running the `run_qc.sh` script, 
so that part of the workflow has already been run and Snakemake will continue 
from that automatically without redoing anything. Once completed you should 
see a bunch of directories and files generated in your current working 
directory, including the `results/` directory containing the final HTML report.
