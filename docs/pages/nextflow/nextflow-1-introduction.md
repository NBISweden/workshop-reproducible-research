<iframe id="iframepdf" src="../../../lectures/nextflow/nextflow.pdf" frameborder="0" width="640" height="480" allowfullscreen="true" mozallowfullscreen="true" webkitallowfullscreen="true"></iframe>

[Nextflow](https://www.nextflow.io/) is a *workflow management system* (WfMS),
and is one of the most common such systems within the bioinformatic and
academic communities. These systems are important for scientific
reproducibility in that they greatly facilitate keeping track of which files
have been processed in what way throughout an entire project.

Nextflow is built from the ground-up to be portable, scalable, reproducible and
usable in a platform-agnostic sense. This means that any workflow you write in
Nextflow can be run locally on your laptop, a computer cluster or a cloud
service (as long as your architecture has the necessary computational
resources). You can also define the compute environment in which each task is
carried out on a per-task basis. You might thus develop your workflow on your
local computer using a minimal test dataset, but run the full analyses with all
samples on *e.g.* a computer cluster. Nextflow can work on both files and
arbitrary values, oftentimes connected in useful and advanced ways.

Nextflow can easily work with dynamic inputs where the exact output is unknown,
*e.g.* the exact number of files or which samples pass some arbitrary quality
control threshold. While Nextflow is based on the Groovy language, you don't
need to know how to code Groovy to be able to write good Nextflow workflows.
Nextflow has a large community centered around it, including the
[nf-core](https://nf-co.re/) curated collection of high quality pipelines used
by *e.g.* the [National Genomics Infrastructure](https://ngisweden.scilifelab.se/).

This tutorial depends on files from the course GitHub repo. Take a look at the
[setup](pre-course-setup) for instructions on how to set it up if you haven't
done so already, then open up a terminal and go to `workshop-reproducible-research/tutorials/nextflow`
and activate your `nextflow-env` Conda environment.
