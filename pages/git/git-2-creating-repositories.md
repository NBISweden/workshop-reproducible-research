In order to create a new git repository, we first need a directory to track.
For this tutorial, go ahead and create a directory called `git_tutorial`, then
navigate into it.

> **Attention!** <br>
> The directory should *not* be within the `workshop-reproducible-research`
> directory, since this is itself a git-tracked directory.

Once we are inside the desired directory, we can *initialise*
git with the following command:

```bash
git init
```

The directory is now a version-tracked directory. How can you know? Run the
command `git status`, which will probably return something like this:

```no-highlight
On branch main

No commits yet

nothing to commit (create/copy files and use "git add" to track)
```

> **Tip** <br>
> If you try to run `git status` in a non-git directory, it will say
> that it is *not a git repository*. The way this works is that git
> adds a hidden directory `.git/` in the root of a git tracked
> directory (run `ls -a` to see it). This hidden directory contains
> all information and settings git needs in order to run and version
> track your files. This also means that your git-tracked directory
> is self-contained, *i.e.* you can simply delete it and everything that
> has to do with git in connection to that directory will be gone.

The text `nothing to commit (create/copy files and use "git add" to track)`
tells us that while we are inside a directory that git is currently tracking,
there are currently no files being tracked; let's add some!

Copy the following files from the `workshop-reproducible-research/git`
directory into your `git_tutorial` directory:

* `Dockerfile`
* `Snakefile`
* `config.yml`
* `environment.yml`

Once you have done that, run `git status` again. It will tell you that there
are files in the directory that are not version tracked by git.

> **Note** <br>
> For the purpose of this tutorial, the exact contents of the files you just
> copied are not important. But you will probably recognize many of them, as
> they are all files used in the MRSA case study described in the
> [introduction to the tutorials](tutorial_intro.md). The details of what
> these files do are described in their respective sessions later in the
> course, but we provide a brief overview here:
> 
> - The `environment.yml` file contains the Conda environment with all the
>   software used in the analysis (see the [Conda tutorial](conda.md)).
> - The `Snakefile` and `config.yml` are both used to define the Snakemake
>   workflow, that we'll go through in the [Snakemake tutorial](snakemake.md).
> - The `Dockerfile` contains the recipe for making a Docker container for
>   the analysis, which will be covered in detail in the
>   [Docker tutorial](docker.md).

> **Quick recap** <br>
> We have used two `git` commands this far:
>
> - `git init` tells git to track the current directory.
> - `git status` is a command you should use *a lot*. It will tell you,
>   amongst other things, the status of your git clone in relation to the
>   online remote repository.
