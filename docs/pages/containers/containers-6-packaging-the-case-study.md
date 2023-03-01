During these tutorials we have been working on a case study about the
multiresistant bacteria MRSA. Here we will build and run a Docker container
that contains all the work we've done so far.

* We've [set up a GitHub repository](git-7-working-remotely) for version control
  and for hosting our project.
* We've defined a [Conda environment](conda-3-projects) that specifies the
  packages we're depending on in the project.
* We've constructed a [Snakemake workflow](snakemake-10-generalizing-workflows)
  that performs the data analysis and keeps track of files and parameters.
* We've written a [R Markdown document](r-markdown-6-the-mrsa-case-study) that
  takes the results from the Snakemake workflow and summarizes them in a report.

The `workshop-reproducible-research/tutorials/containers` directory contains the
final versions of all the files we've generated in the other tutorials:
`environment.yml`, `Snakefile`, `config.yml`, `code/header.tex`, and
`code/supplementary_material.Rmd`. The only difference compared to the other
tutorials is that we have also included the rendering of the Supplementary
Material HTML file into the Snakemake workflow as the rule `make_supplementary`.
Running all of these steps will take some time to execute (around 20 minutes
or so), in particular if you're on a slow internet connection.

Now take a look at `Dockerfile`. Everything should look quite familiar to you,
since it's basically the same steps as in the image we constructed in the
previous section, although some sections have been moved around. The main
difference is that we add the project files needed for executing the workflow
(mentioned in the previous paragraph), and install the conda packages listed in
`environment.yml`. If you look at the `CMD` command you can see that it will
run the whole Snakemake workflow by default.

Now run `docker build` as before, tag the image with `my_docker_project`:

````bash
docker build -t my_docker_project -f Dockerfile .
````

Go get a coffee while the image builds (or
you could use `docker pull nbisweden/workshop-reproducible-research` which
will download the same image).

Validate with `docker image ls`. Now all that remains is to run the whole thing
with `docker run`. We just want to get the results, so mount the directory
`/course/results/` to, say, `mrsa_results` in your current directory.

Well done! You now have an image that allows anyone to exactly reproduce your
analysis workflow (if you first `docker push` to Dockerhub that is).

!!! Tip
    If you've done the [Jupyter tutorial](jupyter-1-introduction), you know that
    Jupyter Notebook runs as a web server. This makes it very well suited for
    running in a Docker container, since we can just expose the port Jupyter
    Notebook uses and redirect it to one of our own. You can then work with the
    notebooks in your browser just as you've done before, while it's actually
    running in the container. This means you could package your data, scripts
    and environment in a Docker image that also runs a Jupyter Notebook server.
    If you make this image available, say on Dockerhub, other researchers could
    then download it and interact with your data/code via the fancy interactive
    Jupyter notebooks that you have prepared for them. We haven't made any
    fancy notebooks for you, but we *have* set up a Jupyter Notebook server.
    Try it out if you want to (replace the image name with your version if
    you've built it yourself):

    ```bash
    docker run -it -p 8888:8888 nbisweden/workshop-reproducible-research \
    jupyter notebook  --ip=0.0.0.0 --allow-root
    ```
