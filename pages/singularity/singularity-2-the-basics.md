In the [Docker tutorial](docker.md) we started by downloading an Ubuntu image.
Let's see how the same thing is achieved with Singularity:

```bash
singularity pull library://ubuntu
```

This pulls an ubuntu image from the [Singularity library](
https://cloud.sylabs.io/library) (somewhat equivalent to Dockerhub). The first
thing you might have noticed is that this command produces a file
`ubuntu_latest.sif` in the current working directory. Singularity, unlike
Docker, stores its images as a single file. Docker on the other hand uses
layers, which can be shared between multiple images, and thus stores downloaded
images centrally (remember the `docker images` command?). A Singularity image
file is self-contained (no shared layers) and can be moved around and shared
like any other file.

To run a command in a Singularity container (equivalent of *e.g.* `docker run
ubuntu uname -a`) we can execute:

```bash
singularity exec ubuntu_latest.sif uname -a
```

This should result in something similar to:
```no-highlight
Linux (none) 4.19.10 #1 SMP Mon Apr 8 00:07:40 CDT 2019 x86_64 x86_64 x86_64 GNU/Linux
```

Now, try to also run the following commands in the ubuntu container in the same
manner as above:

* `whoami`
* `ls -lh`

Notice anything unexpected or different from what you learnt from the Docker
tutorial?

Unlike Docker, Singularity attempts to map parts of your local file system to
the image. By default Singularity bind mounts `$HOME`, `/tmp`, and `$PWD` (the
current working directory) into your container. Also, inside a Singularity
container, you are the same user as you are on the host system.

We can also start an interactive shell (equivalent of *e.g.* `docker run -it
ubuntu`):

```bash
singularity shell ubuntu_latest.sif
```

While running a shell in the container, try executing `pwd` (showing the full
path to your current working directory). See that it appears to be your local
working directory? Try `ls /` to see the files and directory in the root. Exit
the container (`exit`) and run `ls /` to see how it looks on your local system.
Notice the difference?

> **Quick recap** <br>
> In this section we covered:
>
> - How to download a Singularity image using `singularity pull`
> - How to run a command in a Singularity container using `singularity exec`
> - How to start an interactive terminal in a Singularity container using
>   `singularity shell`
