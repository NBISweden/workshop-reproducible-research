Let's assume that you are just about to start a new exciting research project
called *Project A*.

# Creating Conda environments

* Let's make our first Conda environment:

```bash
mamba create -n project_a -c bioconda fastqc
```

This will create an environment called `project_a`, containing FastQC from the
Bioconda channel. Conda will list the packages that will be installed and ask
for your confirmation.

* Once it is done, you can activate the environment:

```bash
mamba activate project_a
```

By default, Conda will add information to your prompt telling you which
environment that is active.

* To see all your environments you can run:

```bash
mamba info --envs
```

The active environment will be marked with an asterisk.

* To see the installed packages and their versions in the active environment,
  run:

```bash
mamba list
```

* To save the installed packages to a file, run:

```bash
mamba env export --from-history > environment.yml
```

Where `--from-history` only reports the packages requested to be installed
and not additional dependencies. A caveat is that if no version was
originally specified, then it is not included in the export file either.

* Now, deactivate the environment by running `mamba deactivate`.
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

* Now, let's add another package (*MultiQC*) to our environment using `conda
  install`. Make sure that `project_a` is the active environment first.

```bash
mamba install -c bioconda multiqc
```

* If we don't specify the package version, the latest available version will be
  installed. What version of MultiQC got installed?
* Run the following to see what versions are available:

```bash
mamba search -c bioconda multiqc
```

* Now try to install a different version of MultiQC, *e.g.*:

```bash
mamba install -c bioconda multiqc=1.13
```

Read the information that Conda displays in the terminal. It probably asks if
you want to downgrade the initial MultiQC installation to the one specified
here (`1.13` in the example). You can only have one version of a given package
in a given environment.

Let's assume that you will have sequencing data in your Project A, and want to
use the latest BBMap software to align your reads.

* Find out what versions of BBMap are available in the Bioconda channel using
  `mamba search -c bioconda bbmap`.
* Now install the *latest* available version of BBMap in your `project_a`
  environment.

Let's further assume that you have an old project (called *Project Old*) where
you know you used BBMap `37.10`. You just got back reviewer comments and they
want you to include some alignment statistics. Unfortunately, you haven't saved
that information so you will have to rerun the alignment. Now, it is essential
that you use the same version of BBMap that your results are based on,
otherwise the alignment statistics will be misleading. Using Conda environments
this becomes simple. You can just have a separate environment for your old
project where you have an old version of BBMap without interfering with your
new Project A where you want the latest version.

* Make a new environment for your old project:

```bash
mamba create -n project_old -c bioconda bbmap=37.10
```

* List your environments (do you remember the command?).
* Activate `project_old` and check the BBMap version (`bbmap.sh --version`).
* Activate `project_a` again and check the BBMap version.

# Removing packages

Now let's try to remove an installed package from the active environment:

```
mamba remove multiqc
```

* Run `mamba deactivate` to exit your active environment.
* Now, let's remove an environment:

```bash
mamba env remove -n project_old
```

After making a few different environments and installing a bunch of packages,
Conda can take up some disk space. You can remove unnecessary files with the
command:

```bash
mamba clean -a
```

This will remove package tar-balls that are left from package installations,
unused packages (*i.e.* those not present in any environments), and cached
data.

> **Quick recap** <br>
> In this section we've learned:
>
> - How to use `mamba install` for installing packages on the fly.
> - How to create, activate and change between environments.
> - How to remove packages or environments and clean up.
