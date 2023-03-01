All of the tutorials and the material in them is dependent on the GitHub
repository for the course. The first step of the setup is thus to download all
the files that you will need, which is done differently depending on which
operating system you have.

At the last day, you will have the opportunity to try out the different tools
on one of your own projects. In case you don't want to use a project you are
currently working on, we have prepared a small-scale project for you. If you
would like to work on your own project, it would be great if you could have the
code and data ready before the end of the course so that you have more time for the exercise.
In case your analysis project contains computationally intense steps it may be
good to scale them down for the sake of the exercise. You might, for example,
subset your raw data to only contain a minuscule part of its original size.

## Setup for Mac / Linux users

First, `cd` into a directory on your computer (or create one) where it makes
sense to download the course directory.

```bash
cd /path/to/your/directory
git clone https://github.com/NBISweden/workshop-reproducible-research.git
cd workshop-reproducible-research
```

!!! tip
    If you want to revisit the material from an older instance of this course,
    you can do that using `git checkout tags/<tag-name>`, e.g.
    `git checkout tags/course_1905`. To list all available tags, use `git tag`.
    Run this command after you have `cd` into `workshop-reproducible-research`
    as described above. If you do that, you probably also want to view the
    same older version of this website. Until spring 2021, the website was
    hosted at https://nbis-reproducible-research.readthedocs.io.
    Locate the version box in the bottom right corner of the website and
    select the corresponding version.

## Setup for Windows users

Using a Windows computer for bioinformatic work has sadly not been ideal most of
the time, but large advanced in recent years have made this quite feasible
through the Windows 10 *Linux subsystem*. This is the only setup for Windows
users that we allow for participants of this course, as all the material has
been created and tested to work on Unix-based systems.

Using the Linux subsystem will give you access to a full command-line bash shell
based on Linux on your Windows 10 PC. For the difference between the Linux Bash
Shell and the PowerShell on Windows 10, see *e.g.* [this article](
https://searchitoperations.techtarget.com/tip/On-Windows-PowerShell-vs-Bash-comparison-gets-interesting).

Install Bash on Windows 10, follow the instructions at *e.g.* one of these
resources:

- [Installing the Windows Subsystem and the Linux Bash](https://docs.microsoft.com/en-us/windows/wsl/install-win10)
- [Installing and using Linux Bash on Windows](https://www.howtogeek.com/249966/how-to-install-and-use-the-linux-bash-shell-on-windows-10/)
- [Installing Linux Bash on Windows](https://itsfoss.com/install-bash-on-windows/)

!!! note

    If you run into error messages when trying to download files through the Linux
    shell (_e.g._ `curl:(6) Could not resolve host`) then try adding the Google
    nameserver to the internet configuration by running `sudo nano /etc/resolv.conf`
    then add `nameserver 8.8.8.8` to the bottom of the file and save it.

Open a bash shell Linux terminal and clone the GitHub repository containing all
files you will need for completing the tutorials as follows. First, `cd` into
a directory on your computer (or create one) where it makes sense to download
the course directory.

!!! tip

    You can find the directory where the Linux distribution is storing all its
    files by typing `explorer.exe .`. This will launch the Windows File Explorer
    showing the current Linux directory.
    Alternatively, you can find the Windows C drive from within the bash shell
    Linux terminal by navigating to `/mnt/c/`.

```bash
cd /path/to/your/directory
git clone https://github.com/NBISweden/workshop-reproducible-research.git
cd workshop-reproducible-research
```

Whenever a setup instruction specifies Mac or Linux (*i.e.* only those two,
with no alternative for Windows), please follow the Linux instructions.

!!! tip

    If you want to revisit the material from an older instance of this course,
    you can do that using `git checkout tags/<tag-name>`, *e.g.* `git checkout
    tags/course_1905`. To list all available tags, use `git tag`. Run this
    command after you have `cd` into `workshop-reproducible-research` as
    described above. If you do that, you probably also want to view the
    same older version of this website. Until spring 2021, the website was
    hosted at https://nbis-reproducible-research.readthedocs.io/en/latest/.
    Locate the version box in the bottom right corner of the website and
    select the corresponding version.

## Installing Git

Chances are that you already have git installed on your computer. You can check
by running *e.g.* `git --version`. If you don't have git, install it following
the instructions [here]( https://git-scm.com/book/en/v2/Getting-Started-Installing-Git).
If you have a very old version of git you might want to update to a later version.

### Configure git

If it is the first time you use git on your computer, you may want to configure
it so that it is aware of your username and email. These should match those that
you have registered on GitHub. This will make it easier when you want to sync
local changes with your remote GitHub repository.

```
git config --global user.name "Mona Lisa"
git config --global user.email "mona_lisa@gmail.com"
```

!!! tip

    If you have several accounts (*e.g.* both a GitHub and Bitbucket account),
    and thereby several different usernames, you can configure git on
    a per-repository level. Change directory into the relevant local git
    repository and run `git config user.name "Mona Lisa"`. This will set the
    default username for that repository only.

You will also need to configure the default branch name to be `main` instead of
`master`:

```bash
git config --global init.defaultBranch "main"
```

The short version of why you need to do this is that GitHub uses `main` as the
default branch while Git itself is still using `master`; please read the box
below for more information.

!!! note

    The default branch name for Git and many of the online resources for hosting
    Git repositories has traditionally been `master`, which historically comes
    from the "master/slave" repositories of
    [BitKeeper](https://mail.gnome.org/archives/desktop-devel-list/2019-May/msg00066.html).
    This has been heavily discussed and in 2020 the decision was made by  many
    ([including GitHub](https://sfconservancy.org/news/2020/jun/23/gitbranchname/))
    to start using `main` instead. Any repository created with GitHub uses this
    new naming scheme since October of 2020, and Git itself is currently
    discussing implementing a similar change. Git did, however, introduce the
    ability to set the default branch name when using `git init` in
    [version 2.28](https://github.blog/2020-07-27-highlights-from-git-2-28/#introducing-init-defaultbranch),
    instead of using a hard-coded `master`. We at NBIS want to be a part of this
    change, so we have chosen to use `main` for this course.

## Installing Conda

Conda is installed by downloading and executing an installer from the Conda
website, but which version you need depends on your operating system:

```bash
# Install Miniconda3 for 64-bit Mac
curl -L https://repo.continuum.io/miniconda/Miniconda3-4.7.12.1-MacOSX-x86_64.sh -O
bash Miniconda3-4.7.12.1-MacOSX-x86_64.sh
rm Miniconda3-4.7.12.1-MacOSX-x86_64.sh
```

```bash
# Install Miniconda3 for 64-bit Linux
curl -L https://repo.continuum.io/miniconda/Miniconda3-4.7.12.1-Linux-x86_64.sh -O
bash Miniconda3-4.7.12.1-Linux-x86_64.sh
rm Miniconda3-4.7.12.1-Linux-x86_64.sh
```

!!! warning

    If you already have installed Conda but want to update, you should be able
    to simply run `conda update conda` and subsequently `conda init`, and skip
    the installation instructions below.

!!! warning
    If you have a newer Apple computer with an M1 chip, make sure you have
    installed [Rosetta](https://support.apple.com/en-us/HT211861) before you run
    the installer. If you want to more fully utilise the new architecture,
    head over to [Miniforge](https://github.com/conda-forge/miniforge#miniforge3)!

The installer will ask you questions during the installation:

- Do you accept the license terms? (Yes)
- Do you accept the installation path or do you want to choose a different one?
  (Probably yes)
- Do you want to run `conda init` to setup Conda on your system? (Yes)

Restart your shell so that the settings in `~/.bashrc`/`~/.bash_profile` can take
effect. You can verify that the installation worked by running:

```bash
conda --version
```

> **Different Condas** <br>
> There are three Conda-related things you may have encountered: the first is
> *Conda*, the package and environment manager we've been talking about so far.
> Second is *Miniconda*, which is the installer for Conda. The third is
> *Anaconda*, which is a distribution of not only Conda, but also over 150
> scientific Python packages. It's generally better to stick with only Conda,
> *i.e.* installing with Miniconda, rather than installing 3 GB worth of
> packages you may not even use.

### Configuring Conda

Lastly, we will setup the default channels (from where packages will be searched
for and downloaded if no channel is specified).

```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

## Installing Snakemake

We will use Conda environments for the set up of this tutorial, but don't worry
if you don't understand exactly what everything does - you'll learn all the
details at the course. First make sure you're currently situated inside the
tutorials directory (`workshop-reproducible-research/tutorials`) and then
create the Conda environment like so:

```bash
conda env create -f snakemake/environment.yml -n snakemake-env
conda activate snakemake-env
```

Check that Snakemake is installed correctly, for example by executing
`snakemake --help`. This should output a list of available Snakemake settings.
If you get `bash: snakemake: command not found` then you need to go back and
ensure that the Conda steps were successful. Once you've successfully completed
the above steps you can deactivate your Conda environment using `conda deactivate`
and continue with the setup for the other tools.

## Installing Nextflow

We'll use Conda to install Nextflow as well: navigate to
`workshop-reproducible-research/tutorials` and create the Conda environment:

```bash
conda env create -f nextflow/environment.yml -n nextflow-env
conda activate nextflow-env
```

Check that Nextflow was installed correctly by running `nextflow -version`. Once
you've successfully completed the installation you can deactive the environment
using `conda deactivate` and continue with the other setups, as needed.

## Installing R Markdown

We also use Conda to install R Markdown: make sure your working directory is in
the tutorials directory (`workshop-reproducible-research/tutorials`) and
install the necessary R packages defined in the `environment.yml`:

```bash
conda env create -f rmarkdown/environment.yml -n rmarkdown-env
```

You can then activate the environment followed by running RStudio in the
background from the command line:

```bash
conda activate rmarkdown-env
rstudio &
```

> **The sluggishness of Conda** <br>
> Some environments are inherently quite complicated in that they have many
> and varied dependencies, meaning that the search space for the entire
> dependency hierarchy becomes huge - leading to slow and sluggish
> installations. This is often the case for R environments. This can be
> improved by using Mamba, a faster wrapper around Conda. Simply run `conda
> install -n base mamba` to install Mamba in your base environment, and
> replace any `conda` command with `mamba` - except activating and
> deactivating environments, which still needs to be done using Conda.

Once you've successfully completed the above steps you can deactivate your Conda
environment using `conda deactivate` and continue with the setup for the other
tools.

> **Windows users** <br>
> In case you are having trouble installing R and RStudio using Conda, both
> run well directly on Windows and you may therefore want to install Windows
> versions of these software for this tutorial (if you haven't done so already).
> Conda is, however, the recommended way. If you're having issues with
> graphical applications, please have a look at [this website](https://seanthegeek.net/234/graphical-linux-applications-bash-ubuntu-windows/);
> scroll down to the "Graphical applications".

> **RStudio and Conda** <br>
> In some cases RStudio doesn't play well with Conda due to differing
> libpaths. The first and simplest thing to try is to always start RStudio from
> the command line (`rstudio &`). If you're still having issues, check the
> available library path by `.libPaths()` to make sure that it points to a path
> within your Conda environment. It might be that `.libPaths()` shows multiple
> library paths, in which case R packages will be searched for by R in all these
> locations. This means that your R session will not be completely isolated in
> your Conda environment and that something that works for you might not work
> for someone else using the same Conda environment, simply because you had
> additional packages installed in the second library location. One way to force
> R to just use the conda library path is to add a `.Renviron` file to the
> directory where you start R with these lines:

    ```
    R_LIBS_USER=""
    R_LIBS=""
    ```

> ... and restart RStudio. The `rmarkdown/` directory in the course materials
> already contains this file, so you shouldn't have to add this yourself, but
> we mention it here for your future projects.

## Installing Jupyter

Let's continue using Conda for installing software, since it's so convenient to
do so! Move in the tutorials directory (`workshop-reproducible-research/tutorials`),
create a Conda environment from the `jupyter/environment.yml` file and test
the installation of Jupyter, like so:

```bash
conda env create -f jupyter/environment.yml -n jupyter-env
conda activate jupyter-env
```

Once you've successfully completed the above steps you can deactivate your Conda
environment using `conda deactivate` and continue with the setup for the other
tools.

## Installing Docker

Installing Docker is quite straightforward on Mac or Windows and a little more
cumbersome on Linux. Note that Docker runs as root, which means that you have to
have `sudo` privileges on your computer in order to install or run Docker. When
you have finished installing docker, regardless of which OS you are on, please
type `docker --version` to verify that the installation was successful!

### macOS

Go to [docker.com](https://docs.docker.com/docker-for-mac/install/#download-docker-for-mac)
and select download option that is suitable for your computer's architecture
(*i.e.* if you have an Intel chip or a newer Apple M1 chip). This will download
a `dmg` file - click on it when it's done to start the installation. This will
open up a window where you can drag the Docker.app to Applications. Close the
window and click the Docker app from the Applications menu. Now it's basically
just to click "next" a couple of times and we should be good to go. You can
find the Docker icon in the menu bar in the upper right part of the screen.

### Linux

How to install Docker differs a bit depending on your Linux distribution, but
the steps are the same. Please follow the instructions for your distribution on
[https://docs.docker.com/engine/install/#server](https://docs.docker.com/engine/install/#server).

!!! Tip
    As mentioned before, Docker needs to run as root. You can achieve this by
    prepending all Docker commands with `sudo`. This is the approach that we
    will take in this tutorial, since the set up becomes a little simpler that way.
    If you plan on continuing using Docker you can get rid of this by adding your
    user to the group `docker`. Here are instructions for how to do this:
    [https://docs.docker.com/engine/installation/linux/linux-postinstall/](https://docs.docker.> com/engine/installation/linux/linux-postinstall/).

### Windows

In order to run Docker on Windows your computer must support *Hardware
Virtualization Technology* and virtualization must be enabled. This is
typically done in BIOS. Setting this is outside the scope of this tutorial,
so we'll simply go ahead as if though it's enabled and hope that it works.

On Windows 10 we will install Docker for Windows, which is available at
[docker.com](https://docs.docker.com/docker-for-windows/install/#download-docker-for-windows).
Click the link *Download from Docker Hub*, and select *Get Docker*. Once the
download is complete, execute the file and follow the
[instructions](https://docs.docker.com/docker-for-windows/install/#install-docker-desktop-on-windows).
You can now start Docker from the Start menu. You can search for it if you
cannot find it; the Docker whale icon should appear in the task bar.

You will probably need to enable integration with the Linux subsystem, if you
haven't done so during the installation of Docker Desktop. Right-click on the
Docker whale icon in the task bar and select *Settings*. Choose *Resources* and
select *WPS integration*. Enable integration with the Linux subsystem and click
*Apply & Restart*; also restart the Linux subsystem.

## Installing Singularity

Installation of Singularity depends, again, on your operating system. When you
have finished, regardless of your OS, please type `singularity --version` to
verify that your installation was successful!

Both Mac and Windows utilise Vagrant, for which the information in the box
below may help you.

> **Vagrant and VirtualBox** <br>
> The Vagrant VirtualBox with Singularity can be started like this:
>
> - Move into the folder `vm-singularity` where you installed Singularity.
> - Type `vagrant up` and once this has finished, verify that the Vagrant
>   VirtualBox is running with `vagrant status`.
> - Now, type `vagrant ssh`, which will open the Vagrant VirtualBox.
> - The first time you open the Vagrant VirtualBox like this, you will have to
>   download the course material to obtain a copy for the Singularity tutorial
>   within the Vagrant VirtualBox by typing `git clone https://github.com/NBISweden/workshop-reproducible-research.git`.

### macOS

Please follow the Mac-specific instructions at the [Singularity website](https://sylabs.io/guides/latest/admin-guide/installation.html#installation-on-windows-or-mac).

### Linux

Follow the Linux-specific instruction at the [Singularity website](https://sylabs.io/guides/latest/admin-guide/installation.html#installation-on-linux).

### Windows

Please follow the Windows-specific instructions at the [Singularity website](https://sylabs.io/guides/latest/admin-guide/installation.html#installation-on-windows-or-mac).

!!! Notes
    Last time we checked, the software "Vagrant Manager" was not available for
    download but the installation of Singularity was successful even without it.

    Version 6.1.28 of "Virtual box for Windows" may not work, please install
    version 6.1.26 from [here](https://www.virtualbox.org/wiki/Download_Old_Builds_6_1)
    in case you encounter problems when trying to start the Vagrant VirtualBox.

## Testing sra-tools
On some computers we've found that the package `sra-tools` which is used in the
course is not working properly. The error seems to be related to some certificate
used to communicate with remote read archives and may affect all environments
with `sra-tools` on the dependency list.

If you run into errors with the program `fastq-dump` from the `sra-tools` package
try the following:

1. Remove `sra-tools` from the relevant environment: `conda remove sra-tools`
2. Download the most recent binaries for your operating system from [here](https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit#the-sra-toolkit-provides-64-bit-binary-installations-for-the-ubuntu-and-centos-linux-distributions-for-mac-os-x-and-for-windows) (example shown for Mac OSX): `curl --output sratoolkit.tar.gz https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-mac64.tar.gz`
3. Create a temporary directory for the installation: `mkdir tmp_out`
4. Extract the binary files: `tar -C tmp_out -zxvf sratoolkit.tar.gz */bin/*`
5. Copy binary files into the conda environment: `cp -r tmp_out/*/bin/* $CONDA_PREFIX/bin/`
6. Remove the downloaded files: `rm -r sratoolkit.tar.gz tmp_out/`
