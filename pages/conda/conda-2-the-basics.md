Let's assume that you are just about to start a new exciting research project
called *Project A*.

# Creating Conda environments

* Let's make our first Conda environment:

```bash
conda create -n project_a -c bioconda fastqc
```

This will create an environment called `project_a`, containing FastQC from the
Bioconda channel. Conda will list the packages that will be installed and ask
for your confirmation.

* Once it is done, you can activate the environment:

```bash
conda activate project_a
```

By default, Conda will add information to your prompt telling you which
environment that is active.

* To see all your environments you can run:

```bash
conda info --envs
```

The active environment will be marked with an asterisk.

* To see the installed packages and their versions in the active environment,
  run:

```bash
conda list
```

* To save the installed packages to a file, run:

```bash
conda env export --from-history > environment.yml
```

where `--from-history` only reports the packages requested to be installed
and not additional dependancies. A caveat is that if no version was
originally specified, then it is not included in the export file either.

* Now, deactivate the environment by running `conda deactivate`.
* List all environments again. Which environment is now marked as active?
* Try to run FastQC:

```bash
fastqc --version
```

* Did it work? Activate your `project_a` environment and run the `fastqc
  --version` command again. Does it work now?

Hopefully the FastQC software was not found in your base environment (unless
you had installed it previously), but worked once your environment was
activated.

# Adding more packages

* Now, let's add another package (*SRA-Tools*) to our environment using `conda
  install`. Make sure that `project_a` is the active environment first.

```bash
conda install -c bioconda sra-tools
```

* If we don't specify the package version, the latest available version will be
  installed. What version of SRA-Tools got installed?
* Run the following to see what versions are available:

```bash
conda search -c bioconda sra-tools
```

* Now try to install a different version of SRA-Tools, *e.g.*:

```bash
conda install -c bioconda sra-tools=2.7.0
```

Read the information that Conda displays in the terminal. It probably asks if
you want to downgrade the initial SRA-Tools installation to the one specified
here (`2.7.0` in the example). You can only have one version of a given package
in a given environment.

Let's assume that you will have sequencing data in your Project A, and want to
use the latest Bowtie2 software to align your reads.

* Find out what versions of Bowtie2 are available in the Bioconda channel using
  `conda search -c bioconda`.
* Now install the *latest* available version of Bowtie2 in your `project_a`
  environment.

Let's further assume that you have an old project (called *Project Old*) where
you know you used Bowtie2 `2.2.5`. You just got back reviewer comments and they
want you to include some alignment statistics. Unfortunately, you haven't saved
that information so you will have to rerun the alignment. Now, it is essential
that you use the same version of Bowtie that your results are based on,
otherwise the alignment statistics will be misleading. Using Conda environments
this becomes simple. You can just have a separate environment for your old
project where you have an old version of Bowtie2 without interfering with your
new Project A where you want the latest version.

* Make a new environment for your old project:

```bash
conda create -n project_old -c bioconda bowtie2=2.2.5
```

* List your environments (do you remember the command?).
* Activate `project_old` and check the Bowtie2 version (`bowtie2 --version`).
* Activate `project_a` again and check the Bowtie2 version.

# Removing packages

Now let's try to remove an installed package from the active environment:

```
conda remove sra-tools
```

* Run `conda deactivate` to exit your active environment.
* Now, let's remove an environment:

```bash
conda env remove -n project_old
```

After making a few different environments and installing a bunch of packages,
Conda can take up some disk space. You can remove unnecessary files with the
command:

```bash
conda clean -a
```

This will remove package tar-balls that are left from package installations,
unused packages (*i.e.* those not present in any environments), and cached
data.

!!! Success "Quick recap"
    In this section we've learned:

    - How to use `conda install` for installing packages on the fly.
    - How to create, activate and change between environments.
    - How to remove packages or environments and clean up.
