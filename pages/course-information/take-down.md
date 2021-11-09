There might be a lot of files stored on your computer after you've taken the
course, depending on how many modules you've gone through. Here are instructions
for how to remove them.

All the tutorials depend on you cloning the `workshop-reproducible-research`
GitHub repo. This can be removed like any other directory; via Finder,
Explorer or `rm -rf workshop-reproducible-research`. Note that this will also
delete the hidden directories `.git`, which contains the history of the repo,
and `.snakemake`, which contains the history of any Snakemake runs.

## Conda

Several of the tutorials use Conda for installing packages. This amounts to
about 2.6 GB if you've done all the tutorials. If you plan on using Conda in
the future you can remove just the packages, or you can remove everything
including Conda itself. 

In order to remove all your Conda environments, you first need to list them:

```bash
conda env list
```

For each of the environments except "base" run the following:

```bash
conda remove -n envname --all
```

And, finally:

```bash
conda clean --all
```

If you also want to remove Conda itself (*i.e.* removing all traces of Conda),
you need to check where Conda is installed. Look for the row "base environment".

```bash
conda info
```

This should say something like `/Users/<user>/miniconda3`. Then remove the
entire Conda directory:

```
rm -rf /Users/<user>/miniconda3
```

Lastly, open your `~/.bashrc` file (or `~/.bash_profile` if on Mac) in a text
editor and remove the path to Conda from PATH.

## Snakemake

Snakemake is installed via Conda and will be removed if you follow the
instructions in the Conda section above. Note that Snakemake also generates
a hidden `.snakemake` directory in the directory where it's run. You can remove
this with the following:

```bash
rm -rf workshop-reproducible-research/tutorials/snakemake/.snakemake
```

## Nextflow

Since we installed Nextflow using Conda we can remove it in the same way as
above. You may also want to remove the `results/` and `work/` directories, which
you can do like so:

```bash
rm -rf workshop-reproducible-research/tutorials/nextflow/results
rm -rf workshop-reproducible-research/tutorials/nextflow/work
```

## Jupyter

Jupyter is installed via Conda and will be removed if you follow the
instructions in the Conda section above.

## Docker

Docker is infamous for quickly taking up huge amounts of space, and some
maintenance is necessary every now and then. Here is how to uninstall Docker
completely. Let's start by removing individual images and containers:

```bash
# Remove unused images
docker image prune

# Remove stopped containers
docker container prune

# Remove unused volumes (not used here, but included for reference)
docker volume prune

# Stop and remove ALL containers
docker container rm $(docker container ls -a -q)

# Remove ALL images
docker image rm $(docker image ls -a -q)
```

Removing Docker itself works differently on the three operating systems, which
is described below:

#### macOS

Click the Docker icon in the menu bar (upper right part of the screen) and
select "Preferences". In the upper right corner, you should find a little bug icon.
Click on that icon and select "Reset to factory defaults". You may have to fill 
in your password. Then select "Uninstall". Once it's done uninstalling, drag the 
Docker app from Applications to Trash.

#### Linux

If you've installed Docker with `apt-get`, uninstall it like this:

```bash
apt-get purge docker-ce
```

Images, containers, and volumes are not automatically removed. To delete all of
them:

```bash
rm -rf /var/lib/docker
```

#### Windows

Uninstall Docker for Windows (on Windows 10) or Docker Toolbox (on Windows 7)
via Control Panel > Programs > Programs and Features. Docker Toolbox will also
have installed Oracle VM VirtualBox, so uninstall that as well if you're not
using it for other purposes.

## Singularity

Singularity images are files that can simply be deleted. Singularity also
creates a hidden directory `.singularity` in your home directory that contains
its cache, which you may delete.

#### Linux

If you want to uninstall Singularity, its removal will depend on the installation 
method you chose. Please refer to the 
[Singularity installation guide](https://sylabs.io/guides/3.4/user-guide/installation.html?highlight=uninstall#install-on-linux) 
to find out which tools and dependencies you need to uninstall.

#### Mac

The DMG file you downloaded to install Singularity Desktop on your Mac contains
an uninstall tool that you can use.

#### Windows

On Windows, you will need to uninstall Git for Windows, VirtualBox, Vagrant and 
Vagrant Manager (see the [Singularity installation guide](https://sylabs.io/guides/3.4/user-guide/installation.html#install-on-windows-or-mac)).
