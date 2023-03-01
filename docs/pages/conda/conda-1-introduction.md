<iframe id="iframepdf" src="../../../lectures/conda/conda.pdf" frameborder="0" width="640" height="480" allowfullscreen="true" mozallowfullscreen="true" webkitallowfullscreen="true"></iframe> 

Conda is a package and environment manager. As a package manager it enables you
to install a wide range of software and tools using one simple command: `conda
install`. As an environment manager it allows you to create and manage multiple
different environments, each with their own set of packages. What are the
benefits of using an environment manager? Some examples include the ability to
easily run different versions of the same package, have different cross-package
dependencies that are otherwise incompatible with each other and, last but not
least, easy installation of all the software needed for an analysis.

Environments are of particular relevance when making bioinformatics projects
reproducible. Full reproducibility requires the ability to recreate the system
that was originally used to generate the results. This can, to a large extent,
be accomplished by using Conda to make a project environment with specific
versions of the packages that are needed in the project. You can read more about
Conda [here](https://conda.io/projects/conda/en/latest/user-guide/concepts/index.html).

A Conda *package* is a compressed tarball (system-level libraries, Python or
other modules, executable programs or other components). Conda keeps track of
the dependencies between packages and platforms - this means that when
installing a given package, all necessary dependencies will also be installed.

Conda packages are typically hosted and downloaded from remote so-called
*channels*. Some widely used channels for general-purpose and bioinformatics
packages are [conda-forge](https://conda-forge.org/) and
[Bioconda](https://bioconda.github.io/), respectively. Both of these are
community-driven projects, so if you're missing some package you can contribute
to the channel by adding the package to it. When installing a Conda package you
specify the package name, version (optional) and channel to download from.

A Conda *environment* is essentially a directory that is added to your PATH and
that contains a specific collection of packages that you have installed.
Packages are symlinked between environments to avoid unnecessary duplication.

This tutorial depends on files from the course GitHub repo. Take a look at the
[setup](pre-course-setup) for instructions on how to set it up if, you haven't
done so already. Then open up a terminal and go to
`workshop-reproducible-research/tutorials/conda`. Instructions below assume
that you are standing in `workshop-reproducible-research/tutorials/conda/`
unless otherwise specified (*e.g.* if it says "create a file", it means save it
in `workshop-reproducible-research/tutorials/conda/`).
