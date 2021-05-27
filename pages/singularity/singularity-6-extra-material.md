# Extra material

A common problem with Singularity is that you can only create local builds if
you are working on a Linux system, as local builds for MacOS and Windows are
currently not supported. This means that you might favour using Docker instead
of Singularity, but what happens when you need to use a HPC cluster such as
Uppmax? Docker won't work there, as it requires root privileges, so Singularity
is the only solution. You can only run Singularity images there, however, not
*build* them...

So, how do you get a Singularity image for use on Uppmax if you can't build it
either locally or on Uppmax? You might think that using remote builds will solve
this, but for a lot of cases this won't help. Since most researchers will want
to work in private Git repositories they can't supply their Conda
`environment.yml` file to remote builds (which only works for public
repositories), which means that you'll have to specify packages manually inside
the container instead.

There is, however, another solution: using Singularity inside Docker. By
creating a bare-bones, Linux-based Docker image with Singularity you can build
Singularity images locally on non-Linux operating systems. This can be either
from Singularity definition files or directly from already existing Docker
images. You can read more about this at the following
[GitHub repository](https://github.com/fasterius/singularity-in-docker).
