# Setup

All of the tutorials and the material in them is dependent on the GitHub
repository for the course. The first step of the setup is thus to download all
the files that you will need, which is done differently depending on which
operating system you have.

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
    same older version of this website. Locate the version box in the bottom 
    right corner of the browser and select the corresponding version.

## Setup for Windows users

There are several different ways to run the course material on a Windows
computer. Neither is perhaps optimal, and the material itself has not been
adapted specifically for Windows. Nevertheless, in principle everything
*should* be possible to run. A few ways you could setup:

- Use the Linux Bash Shell on Windows 10 (see below) _**Recommended for the
  course**_
- Run the course in a Docker container (see below)
- Run as Linux through a virtual machine (and see the Linux setup above)
- Use the Windows 10 PowerShell, install git and clone the course repository
  (see the Mac/Linux setup above)

### Running in the Linux Bash Shell on Windows 10

This will give you access to a full command-line bash shell based on Linux on your
Windows 10 PC. For the difference between the Linux Bash Shell and the PowerShell on Windows
10, see *e.g.* [this article](
https://searchitoperations.techtarget.com/tip/On-Windows-PowerShell-vs-Bash-comparison-gets-interesting).

Install Bash on Windows 10, following the instructions at *e.g.* one of these
resources:

- [Installing the Windows Subsystem and the Linux Bash](https://docs.microsoft.com/en-us/windows/wsl/install-win10)
- [Installing and using Linux Bash on Windows](https://www.howtogeek.com/249966/how-to-install-and-use-the-linux-bash-shell-on-windows-10/)
- [Installing Linux Bash on Windows](https://itsfoss.com/install-bash-on-windows/)

Open a bash shell Linux terminal and clone the GitHub repository
containing all files you will need for completing the tutorials as follows.
First, `cd` into a directory on your computer (or create one) where it makes
sense to download the course directory.

!!! tip
    You can find the directory where the Linux distribution is storing all its files by
    typing `explorer.exe .`. This will launch the Windows File Explorer showing the
    current Linux directory.

```bash
cd /path/to/your/directory
git clone https://github.com/NBISweden/workshop-reproducible-research.git
cd workshop-reproducible-research
```

!!! tip
    If you want to revisit the material from an older instance of this course,
    you can do that using `git checkout tags/<tag-name>`, *e.g.* `git checkout
    tags/course_1905`. To list all available tags, use `git tag`. Run this
    command after you have `cd` into `workshop-reproducible-research` as
    described above. If you do that, you probably also want to view the same
    older version of this website. Locate the version box in the bottom right
    corner of the browser and select the corresponding version.

### Using Docker to run the course

Alternatively, you can use Docker to run the course in a Docker container.
First, open the Windows 10 PowerShell and `cd` into a directory on your computer
(or create one) where it makes sense to download the course directory.
Install Docker by following the instructions in the [Docker tutorial](docker.md#windows).
Then run:

```bash
cd c:/my_dir
docker run -it -p 8888:8888 -v /c/my_dir:/course/ \
    nbisweden/workshop-reproducible-research:slim
```

!!! attention
    Note that we use `/c/my_dir` and not `c:/my_dir` as we normally do on
    Windows. This is required for Docker to parse the command correctly.

This will start an isolated container running Linux, where the directory
`c:/my_dir` is mounted (*i.e.* you can access the files in this Windows
directory within the Linux container, and files edited or created within the
Linux container will appear in this Windows directory). Note that the idea is
that you should edit files in the mounted `c:/my_dir` using an editor in your
normal OS, say Notepad in Windows. The terminal in the container is for running
stuff, not editing.

You should now be at a terminal in the Docker container. Now clone the GitHub
repository containing all the files you will need for the tutorials.

```bash
git clone https://github.com/NBISweden/workshop-reproducible-research.git
cd workshop-reproducible-research
```

#### Setting up persistent Conda environments

Because the root file system of each container is an isolated instance, Conda
environments you create during this course will be lost if you exit the
container or if it is killed for some reason. This means that you will have to
recreate each environment every time you run a new container for
the course. To avoid this, you can make sure that Conda uses a subfolder `envs/`
inside the `/course` directory for storing environments. If you run containers
for the course with some folder on your local machine mounted inside `/course`
that will cause the Conda environments to be available on your local machine
even though they are created inside the container. Should a container be stopped
for some reason you can simply run a new one and activate Conda environments
under `/course/envs`, saving you the trouble of recreating them.

What you have to do is to, after you've started a container with the
`docker run` command above, first create the `envs/` directory:

```bash
mkdir /course/envs
```

Then set the environment variable `CONDA_ENVS_PATH` to `/course/envs`:

```bash
export CONDA_ENVS_PATH="/course/envs"
```

Now when you use Conda to create environments inside the Docker container for
the course, the environments will be placed in `/course/envs` which is also
present on your local system because you mounted the folder inside the
container. **Note that you will have to set the `CONDA_ENVS_PATH` each time you
start a new container for this to work**.

!!! attention
    This usage of persistent Conda environments should be considered an edge
    case of how you use Docker containers. We do this only to make it easier to
    run the course through Docker, but in general we do not advocate creating
    Conda environments separate from the actual Docker container.

Don't worry if you feel that this Docker stuff is a little confusing, it will
become clearer in the [Docker tutorial](docker.md). However, the priority right
now is just to get it running so that you can start working.

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

!!! tip The default branch name
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

!!! attention "Windows users"
    If you are doing the tutorials by running a Docker container on your
    Windows machine, Conda will already be installed for you. You can then jump
    ahead to the last point about setting up the default channels (`conda
    config`) and then go ahead with the practical exercise.

    If you are using the Linux Bash Shell, follow the installation instructions 
    for Linux users.

!!! attention
    If you already have installed Conda but want to update, you should be able
    to simply run `conda update conda` and subsequently `conda init`, and skip
    the installation instructions below.

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

!!! note
    There are three Conda-related things you may have encountered: the first is 
    Conda, the package and environment manager we've been talking about so far. 
    Second is *Miniconda*, which is the installer for Conda. The third is 
    *Anaconda*, which is a distribution of not only Conda, but also over 150 
    scientific Python packages. It's generally better to stick with only Conda, 
    *i.e.* installing with Miniconda, rather than installing 3 GB worth of 
    packages you may not even use.

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
course directory (`workshop-reproducible-research`) and then create the Conda
environment like so:

```bash
conda env create -f snakemake/environment.yml -n snakemake-env
conda activate snakemake-env
```

Check that Snakemake is installed correctly, for example by executing
`snakemake --help`. This should output a list of available Snakemake settings.
If you get `bash: snakemake: command not found` then you need to go back and
ensure that the Conda steps were successful.

!!! attention
    If you look inside `snakemake/environment.yml` you will see that we used the
    package `snakemake-minimal`. This is a slimmed down version that lack some
    features, in particular relating to cloud computing and interacting with
    remote providers such as Google Drive or Dropbox. This was done to speed up
    the installation process. Use the normal `snakemake` package if you need
    those features.

## Installing R Markdown

We also use Conda to install R Markdown: make sure your working directory is in
the course directory (`workshop-reproducible-research`) and install the
necessary R packages defined in the `environment.yml`:

```bash
conda env create -f rmarkdown/environment.yml -n rmarkdown-env
```

You can then activate the environment followed by running RStudio in the
background from the command line:

```bash
conda activate rmarkdown-env
rstudio &
```

!!! note "RStudio and Conda"
    In some cases RStudio doesn't play well with Conda due to differing
    libpaths. To fix this, first check the available library path by
    `.libPaths()` to make sure that it points to a path within your conda
    environment. It might be that `.libPaths()` shows multiple library paths, in
    which case R packages will be searched for by R in all these locations.
    This means that your R session will not be completely isolated in your Conda
    environment and that something that works for you might not work for
    someone else using the same Conda environment, simply because you had
    additional packages installed in the second library location. One way to
    force R to just use the conda library path is to add a `.Renviron` file to
    the directory where you start R with these lines:

    ```
    R_LIBS_USER=""
    R_LIBS=""
    ```

    ... and restart RStudio. The `rmarkdown/` directory in the course materials
    already contains this file, so you shouldn't have to add this yourself, but
    we mention it here for your future projects.

!!! attention "Windows users"
    Although most of the tutorials are best to run in the Linux Bash Shell or
    in a Docker container if you are a Windows user (see information in the
    [intro](tutorial_intro.md)), both R and RStudio run well directly on
    Windows. You may therefore want to install Windows versions of these
    software (if you haven't done so already) when doing this tutorial, if
    you're having trouble using Conda. Conda is, however, the recommended way.

## Installing Jupyter

Let's continue using Conda for installing software, since it's so convenient to
do so! Create an environment from the `jupyter/environment.yml` file and test
the installation of Jupyter, like so:
 
```bash
conda env create -f jupyter/environment.yml -n jupyter-env
conda activate jupyter-env
```

!!! attention "Windows users"
    If you are doing these exercises through a Docker container you also need
    the run the following:
    
    ```bash
    mkdir -p -m 700 /root/.jupyter/ && \
    echo "c.NotebookApp.ip = '0.0.0.0'" >> \
        /root/.jupyter/jupyter_notebook_config.py
    ```

## Installing Docker

Installing Docker is quite straightforward on Mac or Windows and a little more
cumbersome on Linux. Note that Docker runs as root, which means that you have to
have `sudo` privileges on your computer in order to install or run Docker.

### macOS

Go to [docker.com]( https://docs.docker.com/docker-for-mac/install/#download-docker-for-mac)
and select "Get Docker for Mac (Stable)". This will download a `dmg` file. Click
on it once it's done to start the installation. This will open up a window
where you can drag the Docker.app to Applications. Close the window and click
the Docker app from the Applications menu. Now it's basically just to click
"next" a couple of times and we should be good to go. You can find the Docker
icon in the menu bar in the upper right part of the screen.

### Windows

The instructions are different depending on if you have Windows 10 or Windows 7.
In order to run Docker on Windows your computer must support Hardware 
Virtualization Technology and virtualization must be enabled. This is 
typically done in BIOS. Setting this is outside the scope of this tutorial, 
so we'll simply go ahead as if though it's enabled and hope that it works.

On Windows 10 we will install Docker for Windows, which is available at
[docker.com](https://docs.docker.com/docker-for-windows/install/#download-docker-for-windows).
Select "Get Docker for Windows (Stable)".

1. Once it's downloaded, double-click `Docker for Windows Installer.exe` to run
   the installer.

2. Follow the install wizard and accept the license, authorize the installer,
   and proceed with the install. You will be asked to authorize Docker.app with
   your system password during the install process. Click Finish to exit the
   installer.

3. Start Docker from the Start menu. You can search for it if you cannot find
   it. The Docker whale icon should appear in the task bar. Open the Windows 10
   PowerShell to run the tutorial.

4. If you would like to use the Linux Bash Shell instead of the Windows 10
   PowerShell, you might need to enable integration with the Linux app you
   installed. Right-click on the Docker whale icon in the task bar and select
   Settings. Choose Resources and in there, select WPS integration. Enable
   integration with the Linux app you installed and click Apply & Restart.
   Restart also the Linux app.

On Windows 7 we will instead use Docker Toolbox, which is available at
[docker.com](https://docs.docker.com/toolbox/toolbox_install_windows/). Select
"Get Docker Toolbox for Windows".

1. Install Docker Toolbox by double-clicking the installer. Step through the
   installation and accept all the defaults. If Windows security dialog prompts
   you to allow the program to make a change, choose Yes. If you get a prompt
   asking "Would you like to install this device software?" select Install.

2. You should now have a Docker Quickstart icon on the desktop.

### Linux

How to install Docker differs a bit depending on your Linux distribution, but
the steps are the same. For details on how to do it on your distribution see
[https://docs.docker.com/engine/installation/#server](https://docs.docker.com/engine/install/#server).

Here we show how to do it for Ubuntu, which is the most common desktop
distribution. Docker requires a 64-bit Ubuntu version 14.04 or higher. If your
OS is from 2015 or earlier you can double check this with `lsb_release -a`. If
it's newer you're probably fine. The same instructions apply to other 
distributions based on Ubuntu, such as Elementary OS or Linux Mint, but you 
would have to map to the corresponding Ubuntu version and use that instead of 
`$(lsb_release -cs)` below (see 
[here](https://en.wikipedia.org/wiki/Linux_Mint_version_history#Release_history)
for Mint).

1. Add the GPG key for the official Docker repository to the system:

```bash
curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo apt-key add -
```

2. Add the Docker repository to APT sources:

```bash
sudo add-apt-repository "deb [arch=amd64] \
    https://download.docker.com/linux/ubuntu $(lsb_release -cs) stable"
```

3. Update the package database with the Docker packages from the new repo:

```bash
sudo apt-get update
```

4. Install Docker Community Edition:

```bash
sudo apt-get install -y docker-ce
```

5. Docker should now be installed, the daemon started, and the process enabled
   to start on boot. Check that it's running:

```bash
sudo systemctl status docker
```

The output should say something about "Active: active (running) since..".

!!! tip
    As mentioned before, Docker needs to run as root. You can achieve this by
    prepending all Docker commands with `sudo`. This is the approach that we
    will take in this tutorial, since the set up becomes a little simpler that way. 
    If     you plan on continuing using Docker you can get rid of this by adding your
    user to the group `docker`. Here are instructions for how to do this:
    [https://docs.docker.com/engine/installation/linux/linux-postinstall/][].


## Installing Singularity

Installation of Singularity depends, again, on your operating system. Here are
instructions for each of them:

### macOS

Download the Singularity Desktop DMG file from [here](https://sylabs.io/singularity-desktop-macos/)
and follow the instructions. Note that this is a beta version and not all
features are available yet.

!!! attention
    Make sure you that 'Singularity networking' is checked during installation

### Linux

Follow the instructions [here](https://sylabs.io/guides/3.4/user-guide/installation.html#distribution-packages-of-singularity).

### Windows

Installing on Windows requires running Singularity through a Vagrant Box, which
may be tricky. See [instructions here](https://sylabs.io/guides/3.4/user-guide/installation.html#install-on-windows-or-mac).

The Vagrant VirtualBox with Singularity can be started on your Windows 10 PC
like this:

* Open the Vagrant Manager app
* Open the Git Bash and move with `cd` into the folder `vm-singularity` where
  you installed Singularity
* Type `vagrant up` and once this has finished, verify that the Vagrant
  VirtualBox is running with `vagrant status`
* Now, type `vagrant ssh`, which will open the Vagrant VirtualBox
* Here, you will have to download the course material again to obtain a copy
  within the Vagrant VirtualBox for the Singularity tutorial by typing 
  `git clone https://github.com/NBISweden/workshop-reproducible-research.git`
