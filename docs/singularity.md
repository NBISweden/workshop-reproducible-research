# Introduction to Singularity

!!! attention
    This is a beta version tutorial that is currently only aimed at getting a Docker image to run using singularity. Although aiming to be general, it is also somewhat focused on getting this running on the Swedish university compute cluster (HPC) Uppmax.

## Tell me more
* [Singularity homepage](http://singularity.lbl.gov/)
* [Singularity publication in PLoS ONE](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0177459)
* [Singularity hub](https://www.singularity-hub.org) - like Dockerhub for Singularity images
* [Uppmax info](https://www.uppmax.uu.se/support/user-guides/singularity-user-guide/)

# Set up

## Install singularity

Refer to [this documentation](http://singularity.lbl.gov/docs-installation).

!!! note
    Singularity is already installed on Uppmax. For this tutorial, you do not need to install it on your local system.

# Practical exercise

## A quick and dirty intro

Before we start, make a new directory and `cd` into it!  

Now, let's start like we did in the Docker tutorial:

```
singularity pull docker://ubuntu:latest
```

This creates an image file in your working directory (unlike Docker). Note the message *Cache folder set to ...* Keep an eye on this, since it can grow to several Gb fast. You need to clean it manually.

!!! tip
    Alternatively you can set the cache directory to e.g. `/scratch` (or something else) using `export SINGULARITY_CACHEDIR=/scratch`.

Now, (similar to `docker run -it ubuntu`), run:
```
singularity shell ubuntu-latest.img
```

You are now in the container, run e.g.:

```
uname -a
ls
pwd
ls $HOME
touch hello
```

This does not behave as Docker as you can see. It appears as if you are still on the host system. This is because, by default, singularity bind mounts `/home/$USER`, `/tmp`, and `$PWD` into your container at runtime. Additionally, the system administrator may have setup other bind mounts to be used by default (try e.g. `ls /` and see if you recognize any directories).

Nevertheless, you are in the container, run `cat /etc/os-release` and compare the output to that on your host system. Exit the container using `exit`. Notice that the file `hello` exists outside the container (since the working directory was mounted).

In order to actual contain our analysis as much as possible we can use the following:

```
singularity shell --cleanenv --containall ubuntu-latest.img
```

Now, try `pwd` and `ls`. Notice that unlike before you are not in the host working directory but in the container `$HOME` and that directory is empty (i.e. not mounted). Also run `ls /` and compare to what you got previously.

Notice however that this is not completely contained. Your user is kept (`whoami`). This is also reflected in your home directory path in the container (`pwd $HOME`).

You could make a file (`touch hej`) but you would not be able to access this outside the container. Exit the container and verify this.

Even though we normally want to isolate the container from the host, it is useful to bind mount a specific directory. Make a new directory within your working directory (`mkdir results`). Then run:

```
singularity shell --cleanenv --containall --bind results:$HOME/output ubuntu-latest.img
```

This means that `$HOME/output` in the container (`output/` was created when running `singularity shell`) will map to `$PWD/results/` on the host system.

Try adding a file (`touch output/table1.txt`), exit the container, and check the contents of `results/`.

We could also replace the home directory path in the container and mount a local directory of choice. Using `--bind` for this (`--bind results:/home`) will generate a warning and move us to `/` though, since singularity tries to set the container working directory to `$HOME` (i.e. `/home/$(whoami)`) and will not find that anymore. Instead we can use the `--home` option:

```
singularity shell --cleanenv --containall --home results:/home ubuntu-latest.img
```

Run `pwd`, notice we now start in `/home`. Run `ls` and se the `table1.txt` that was already in `results/` on the host.

!!! note
    There are additional commands that we have not covered here. See `singularity --help` for a list.

## Running the MRSA workflow

Get the same image from Dockerhub that we used in the Docker tutorial (you probably want to do this on a compute node rather than the login node if you are working on an HPC cluster):
```
singularity pull -n mrsa_workflow.img docker://scilifelablts/reproducible_research_course
```

Next, shell into the container:

```
singularity shell --cleanenv --containall --bind results/:/dev/shm/results mrsa_workflow.img
```

!!! note
    Here we use a trick and mount the host directory `results/` to `/dev/shm` in the container. This is where we will run the analysis and the reason is that you will otherwise probably run into disk space issues, for some reason. Try it if you want.

The files that we need are in `/files` in the container. Unfortunately this directory, as most, are read-only. Therefore we need to move them before we run the analysis. It would make sense to move them to `$HOME` or similar, however, as mentioned in the note above, this will cause disk space issues while running the analysis. Instead we use the shared memory:

```
cd /dev/shm
cp -r /files/* .
ls
```

!!! tip
    We could also have e.g. git cloned the code on the host and bind mounted that directory into the container, instead of copying files within the container as we do here. The choice is yours!

Now that the files are in place, we can execute the workflow:

```
snakemake --configfile config.yml
ls
```
Notice that a few directories have been created (as you will recognize from the Snakemake tutorial), and most important, the `results/` directory is now filled with result files! Good thing we mounted that directory.

Before we exit, just a side-note. Run `ls $HOME`. See that a directory `ncbi` is created by the workflow. This is a side-effect of the `fastq-dump` program. Now that we ran with `--containall` it is no big deal. However, without that option, the `ncbi` directory would have appeared in you `$HOME` on the host system. Another good reason to always attempt to contain the container as much as possible! Ok, go ahead and exit the container. Check the host `results` directory, where you will find the result files (like `supplementary_material.pdf`) safely stored!

That's all folks.


## Cleaning up

Remember to delete the cache! Defaults to `$HOME/.singularity/`. Unlike Docker, containers are not stored. What remains are the image files that you pulled (`ubuntu-latest.img` and `mrsa_workflow.img`).
