All of the tutorials and the material in them is dependent on the GitHub
repository for the course. The first step of the setup is thus to download all
the files that you will need, which is done differently depending on which
operating system you have.

On the last day, in the _Putting the pieces together_ session we will give
examples of how we use these tools in our day-to-day work. During the course,
spend some time thinking about how these tools could be useful for you in your
own project(s). After you've gone through the tutorial you may feel that some of
the tools could make your life easier from the get-go, while others may take
some time to implement efficiently (and some you may never use again after the
course). Our idea with the _Putting the pieces together_ session is to have an
open discussion about where to go from here.

## Setup for Mac / Linux users

First, `cd` into a directory on your computer (or create one) where it makes
sense to download the course directory.

```bash
cd /path/to/your/directory
git clone https://github.com/NBISweden/workshop-reproducible-research.git
cd workshop-reproducible-research
```

> **Tip** <br>
> If you want to revisit the material from an older instance of this course,
> you can do that using `git switch -d tags/<tag-name>`, e.g.
> `git switch -d tags/course_1905`. To list all available tags, use `git tag`.
> Run this command after you have `cd` into `workshop-reproducible-research`
> as described above. If you do that, you probably also want to view the
> same older version of this website. Until spring 2021, the website was
> hosted at [ReadTheDocs](https://nbis-reproducible-research.readthedocs.io).
> Locate the version box in the bottom right corner of the website and
> select the corresponding version.

## Setup for Windows users

Using a Windows computer for bioinformatic work has sadly not been ideal most of
the time, but large advanced in recent years have made this quite feasible
through the Windows 10 _Linux subsystem_. This is the only setup for Windows
users that we allow for participants of this course, as all the material has
been created and tested to work on Unix-based systems.

Using the Linux subsystem will give you access to a full command-line bash shell
based on Linux on your Windows 10 PC. For the difference between the Linux Bash
Shell and the PowerShell on Windows 10, see _e.g._ [this article](https://searchitoperations.techtarget.com/tip/On-Windows-PowerShell-vs-Bash-comparison-gets-interesting).

Install Bash on Windows 10, follow the instructions at _e.g._ one of these
resources:

- [Installing the Windows Subsystem and the Linux Bash](https://docs.microsoft.com/en-us/windows/wsl/install-win10)
- [Installing and using Linux Bash on Windows](https://www.howtogeek.com/249966/how-to-install-and-use-the-linux-bash-shell-on-windows-10/)
- [Installing Linux Bash on Windows](https://itsfoss.com/install-bash-on-windows/)

> **Note** <br>
> If you run into error messages when trying to download files through the Linux
> shell (_e.g._ `curl:(6) Could not resolve host`) then try adding the Google
> name server to the internet configuration by running `sudo nano
/etc/resolv.conf` then add `nameserver 8.8.8.8` to the bottom of the file and
> save it.

> **Important!** <br>
> Whenever a setup instruction specifies Mac or Linux (_i.e._ only those two, with no alternative for Windows),
> **please follow the Linux instructions.**

Open a bash shell Linux terminal and clone the GitHub repository containing all
files you will need for completing the tutorials as follows. First, `cd` into
a directory on your computer (or create one) where it makes sense to download
the course directory.

> **Tip** <br>
> You can find the directory where the Linux distribution is storing all its
> files by typing `explorer.exe .`. This will launch the Windows File Explorer
> showing the current Linux directory.
> Alternatively, you can find the Windows C drive from within the bash shell
> Linux terminal by navigating to `/mnt/c/`.

```bash
cd /path/to/your/directory
git clone https://github.com/NBISweden/workshop-reproducible-research.git
cd workshop-reproducible-research
```

## Installing Git

Chances are that you already have git installed on your computer. You can check
by running _e.g._ `git --version`. If you don't have git, install it following
the instructions [here](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git).
If you have a very old version of git you might want to update to a later
version. If you're on a Mac you can also install it using [Homebrew](https://brew.sh/)
and simple `brew install git`.

### Configure git

If it is the first time you use git on your computer, you may want to configure
it so that it is aware of your username and email. These should match those that
you have registered on GitHub. This will make it easier when you want to sync
local changes with your remote GitHub repository.

```
git config --global user.name "Mona Lisa"
git config --global user.email "mona_lisa@gmail.com"
```

> **Tip** <br>
> If you have several accounts (_e.g._ both a GitHub and Bitbucket account),
> and thereby several different usernames, you can configure git on
> a per-repository level. Change directory into the relevant local git
> repository and run `git config user.name "Mona Lisa"`. This will set the
> default username for that repository only.

You will also need to configure the default branch name to be `main` instead of
`master`:

```bash
git config --global init.defaultBranch "main"
```

The short version of why you need to do this is that GitHub uses `main` as the
default branch while Git itself is still using `master`; please read the box
below for more information.

> **The default branch name** <br>
> The default branch name for Git and many of the online resources for hosting
> Git repositories has traditionally been `master`, which historically comes
> from the "master/slave" repositories of
> [BitKeeper](https://mail.gnome.org/archives/desktop-devel-list/2019-May/msg00066.html).
> This has been heavily discussed and in 2020 the decision was made by many
> ([including GitHub](https://sfconservancy.org/news/2020/jun/23/gitbranchname/))
> to start using `main` instead. Any repository created with GitHub uses this
> new naming scheme since October of 2020, and Git itself is currently
> discussing implementing a similar change. Git did, however, introduce the
> ability to set the default branch name when using `git init` in
> [version 2.28](https://github.blog/2020-07-27-highlights-from-git-2-28/#introducing-init-defaultbranch),
> instead of using a hard-coded `master`. We at NBIS want to be a part of this
> change, so we have chosen to use `main` for this course.

### GitHub setup

[GitHub](https://github.com) is one of several online hosting platforms for Git
repositories. We'll go through the details regarding how Git and GitHub are
connected in the course itself, so for now we'll stick to setting up your
account and credentials.

If you have not done so already, go to [github.com](https://github.com/join) and
create an account. You can also create an account on another online hosting
service for version control, _e.g._ [Bitbucket](https://bitbucket.org) or
[GitLab](https://about.gitlab.com/). The exercises in this course are written
with examples from GitHub (as that is the most popular platform with the most
extensive features), but the same thing can be done on alternative services,
although the exact menu structure and link placements differ.

Any upload to and from GitHub requires you to authenticate yourself. GitHub
used to allow authentication with your account and password, but this is no
longer the case - using SSH keys is required instead. Knowing exactly what these
are is not necessary to get them working, but we encourage you to read the box
below to learn more about them! GitHub has excellent, platform-specific
instructions both on how to [generate](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent)
and [add](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account)
SSH keys to your account, so please follow those instructions.

> **SSH keys and authentication** <br>
> Using SSH (Secure Shell) for authentication basically entails setting up a
> pair of keys: one private and one public. You keep the private key on your
> local computer and give the public key to anywhere you want to be able to
> connect to, _e.g._ GitHub. The public key can be used to encrypt messages that
> _only_ the corresponding private key can decrypt. A simplified description of
> how SSH authentication works goes like this:
>
> 1. The client (_i.e._ the local computer) sends the ID of the SSH key pair it
>    would like to use for authentication to the server (_e.g._ GitHub)
> 2. If that ID is found, the server generates a random number and encrypts this
>    with the public key and sends it back to the client
> 3. The client decrypts the random number with the private key and sends it
>    back to the server
>
> Notice that the private key always remains on the client's side and is never
> transferred over the connection; the ability to decrypt messages encrypted
> with the public key is enough to ascertain the client's authenticity. This is
> in contrast with using passwords, which are themselves sent across a
> connection (albeit encrypted). It is also important to note that even though
> the keys come in pairs it is impossible to derive the private key from the
> public key. If you want to read more details about how SSH authentication work
> you can check out [this website](https://www.digitalocean.com/community/tutorials/understanding-the-ssh-encryption-and-connection-process),
> which has more in-depth information than we provide here.

## Installing Mamba

### Mamba or Conda?

Maybe you've worked with the Conda package manager before, and you're wondering
what Mamba is? Mamba is, simply put, a faster implementation of Conda. Mamba has
quickly grown and matured to the point that we are almost exclusively using it
in our own daily work rather than Conda - we are thus reflecting this
wide-spread adopting in the course material as well. Conveniently there is
almost no difference in the way the two programs work on the command line. You
simply change `conda` to `mamba` and keep working as you've done before (see the
exception under _Configuring Conda_ below). This also means that if you already
have conda installed you can keep using it for this course, however we strongly
recommend you to try out mamba in order to make your environment managing more
efficient.

You will notice that we still use the terms _Conda environment_, _Conda
packages_ _etc._ throughout the course and that the tutorial pages still have
_Conda_ in the title. This is because Conda is the original package manager with
a widely adopted terminology, and Mamba is a re-implementation of Conda; we hope
this is not too confusing.

### If you already have Mamba installed

If you already have installed Mamba you can make sure you're using the latest
version by running `mamba update -n base mamba` and skip the installation
instructions below.

### If you have Conda installed

conda install mamba -n base mamba -c conda-forge

```

Check that installation worked by running:

```

mamba --version

````

### If you have neither Mamba nor Conda installed

Mamba is installed by downloading and executing a [Mambaforge](https://github.com/conda-forge/miniforge#mambaforge)
installer for your operating system.

```bash
# Install Mambaforge for 64-bit Mac
curl -L https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-MacOSX-x86_64.sh -O
bash Mambaforge-MacOSX-x86_64.sh
rm Mambaforge-MacOSX-x86_64.sh
````

```bash
# Install Mambaforge 3 for 64-bit Mac (Apple chip)
curl -L  https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-MacOSX-arm64.sh -O
bash Mambaforge-MacOSX-arm64.sh
rm Mambaforge-MacOSX-arm64.sh
```

```bash
# Install Mambaforge3 for 64-bit Linux
curl -L https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh -O
bash Mambaforge-Linux-x86_64.sh
rm Mambaforge-Linux-x86_64.sh
```

The installer will ask you questions during the installation:

- Do you accept the license terms? (Yes)
- Do you accept the installation path or do you want to choose a different one?
  (Probably yes)
- Do you want the installer to initialize Mambaforge (Yes)

Restart your shell so that the settings in `~/.bashrc` or `~/.bash_profile` can
take effect. You can verify that the installation worked by running:

```bash
mamba --version
```

> **Important!** <br>
> The Mamba docs specify a couple of things to keep in mind when using Mamba.
> First of all, `mamba` should be installed in the `base` environment and no
> other packages should be installed into `base`. Furthermore, mixing of the
> `conda-forge` and `defaults` channels should be avoided as the default
> Anaconda channels are incompatible with `conda-forge`.

> **Different Mamba/Conda flavours** <br>
> You may come across several flavours of both Mamba and Conda. For Mamba
> there's the _Miniforge_ installer which allows you to install the `mamba`
> command line tool that works as a replacement for `conda`. There's also
> `micromamba`, a small stand-alone C++ program developed mainly for continuous
> integration pipelines. For Conda there's _Miniconda_, which is the installer
> for Conda. The third is _Anaconda_, which is a distribution of not only Conda,
> but also over 150 scientific Python packages. If you want to use Conda it's
> generally better to stick with the Miniconda installation, rather than
> installing 3 GB worth of packages you may not even use.

### Configuring Mamba/Conda

At the moment the `config` subcommand is not implemented in `mamba`. This means
that when you want to configure your mamba or conda installation you still need
to rely on `conda`.

As a last step we will set up the default channels (from where packages will be
searched for and downloaded if no channel is specified):

```
conda config --add channels bioconda
conda config --add channels conda-forge
```

And we will also set so called 'strict' channel priority, which ensures higher
stability and better performance (see details about this setting by running the
following:

```
conda config --set channel_priority strict
```

### Mamba/Conda on new Macs

If you have one of the newer Macs with Apple chips (the M-series) you may run
into some problems with certain Conda packages that have not yet been built for
the ARM64 architecture. The [Rosetta](https://support.apple.com/en-us/HT211861)
software allows ARM64 Macs to use software built for the old AMD64 architecture,
which means you can always fall back on creating AMD/Intel-based environments
and use them in conjunction with Rosetta. This is how you do it:

```bash
CONDA_SUBDIR=osx-64 <mamba-command>
conda activate <env>
conda config --env --set subdir osx-64
```

The first command creates the Intel-based environment, while the last one
makes sure that subsequent commands are also using the Intel architecture. If
you don't want to remember and do this manually each time you want to use
AMD64/Rosetta you can check out [this bash script](https://github.com/fasterius/dotfiles/blob/main/scripts/intel-conda-env.sh).

## Installing Snakemake

We will use Conda environments for the set up of this tutorial, but don't worry
if you don't understand exactly what everything does - you'll learn all the
details at the course. First make sure you're currently situated inside the
tutorials directory (`workshop-reproducible-research/tutorials`) and then create
the Conda environment like so:

```bash
mamba env create -f snakemake/environment.yml -n snakemake-env
mamba activate snakemake-env
```

> **ARM64 users:** <br>
> Some of the packages in this environment is not available for the ARM64
> architecture, so you'll have to follow the [instructions above](#mamba/Conda-on-new-macs).

Check that Snakemake is installed correctly, for example by executing `snakemake
--help`. This should output a list of available Snakemake settings. If you get
`bash: snakemake: command not found` then you need to go back and ensure that
the Mamba steps were successful. Once you've successfully completed the above
steps you can deactivate the environment using `mamba deactivate` and continue
with the setup for the other tools.

## Installing Nextflow

The easiest way to install Nextflow is the official one, which is to just run the
following code:

```bash
curl -s https://get.nextflow.io | bash
```

This will give you the `nextflow` file in your current directly - move this file
to a directory in your `PATH`, _e.g._ `/usr/bin/`.

If you're getting Java-related errors, you can either try to [update your Java
installation](https://www.nextflow.io/docs/latest/getstarted.html#requirements)
(Nextflow requires Java 11 or later) or install Nextflow using Mamba. If you
want to use Mamba, navigate to `workshop-reproducible-research/tutorials` and
create the environment:

```bash
mamba env create -f nextflow/environment.yml -n nextflow-env
mamba activate nextflow-env
```

> **ARM64 users:** <br>
> Some of the packages in this environment is not available for the ARM64
> architecture, so you'll have to follow the [instructions above](#mamba/Conda-on-new-macs).

Check that Nextflow was installed correctly by running `nextflow -version`. If
you successfully installed Nextflow using Mamba you can now deactivate the
environment using `mamba deactivate` and continue with the other setups, as
needed.

## Installing Quarto

Installing Quarto is easiest by going to the [official
website](https://quarto.org/docs/get-started/) and downloading the
OS-appropriate package and following the installation instructions. You also
need to install a LaTeX distribution to be able to render Quarto documents to
PDF, which can be done using Quarto itself:

```bash
quarto install tinytex
```

While we're not installing Quarto _itself_ using Mamba, we _will_ install some
software packages that are used in the Quarto tutorial using Mamba: make sure
your working directory is in the tutorials directory (`workshop-reproducible-research/tutorials`)
and install the necessary packages defined in the `environment.yml`:

```bash
mamba env create -f quarto/environment.yml -n quarto-env
```

> **ARM64 users:** <br>
> Some of the packages in this environment is not available for the ARM64
> architecture, so you'll have to follow the [instructions above](#mamba/Conda-on-new-macs).

## Installing Jupyter

Let's continue using Mamba for installing software, since it's so convenient to
do so! Move into the tutorials directory (`workshop-reproducible-research/tutorials`),
create a Conda environment from the `jupyter/environment.yml` file and test
the installation of Jupyter, like so:

```bash
mamba env create -f jupyter/environment.yml -n jupyter-env
mamba activate jupyter-env
```

Once you've successfully completed the above steps you can deactivate the
environment using `mamba deactivate` and continue with the setup for the other
tools.

## Installing Docker

Installing Docker (specifically Docker Desktop) is quite straightforward on Mac,
Windows and Linux distributions. Note that Docker runs as root, which means that
you have to have `sudo` privileges on your computer in order to install or run
Docker. When you have finished installing docker, regardless of which OS you are
on, please type `docker --version` to verify that the installation was
successful.

> **Docker for older versions of OSX/Windows** <br>
> The latest version of Docker may not work if you have an old version of either
> OSX or Windows. You can find older Docker versions that may be compatible for
> you if you go to
> [https://docs.docker.com/desktop/](https://docs.docker.com/desktop/) and click
> "Previous versions" in the left side menu.

### MacOS

Go to [docker.com](https://docs.docker.com/docker-for-mac/install/#download-docker-for-mac)
and select the download option that is suitable for your computer's architecture
(_i.e._ if you have an Intel chip or a newer Apple M1 chip). This will download
a `dmg` file - click on it when it's done to start the installation. This will
open up a window where you can drag the Docker.app to Applications. Close the
window and click the Docker app from the Applications menu. Now it's basically
just to click "next" a couple of times and we should be good to go. You can find
the Docker icon in the menu bar in the upper right part of the screen.

### Linux

Go to the [linux-install](https://docs.docker.com/desktop/install/linux-install/)
section of the Docker documentation and make sure that your computer meets the
system requirements. There you can also find instructions for different Linux
distributions in the left sidebar under _Installation per Linux distro_.

### Windows

In order to run Docker on Windows your computer must support _Hardware
Virtualisation Technology_ and virtualisation must be enabled. This is typically
done in BIOS. Setting this is outside the scope of this tutorial, so we'll
simply go ahead as if though it's enabled and hope that it works.

On Windows 10 we will install Docker for Windows, which is available at
[docker.com](https://docs.docker.com/docker-for-windows/install/#download-docker-for-windows).
Click the link _Download from Docker Hub_, and select _Get Docker_. Once the
download is complete, execute the file and follow the [instructions](https://docs.docker.com/docker-for-windows/install/#install-docker-desktop-on-windows).
You can now start Docker from the Start menu. You can search for it if you
cannot find it; the Docker whale icon should appear in the task bar.

You will probably need to enable integration with the Linux subsystem, if you
haven't done so during the installation of Docker Desktop. Right-click on the
Docker whale icon in the task bar and select _Settings_. Choose _Resources_ and
select _WPS integration_. Enable integration with the Linux subsystem and click
_Apply & Restart_; also restart the Linux subsystem.
