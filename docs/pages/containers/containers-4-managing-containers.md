When you start a container with `docker run` it is given an unique id that you
can use for interacting with the container. Let's try to run a container from
the image we just created:

```bash
docker run my_docker_conda
```

If everything worked `run_qc.sh` is executed and will first download and then
analyse the three samples. Once it's finished you can list all containers,
including those that have exited.

```bash
docker container ls --all
```

This should show information about the container that we just ran. Similar to:

```
CONTAINER ID        IMAGE                    COMMAND                  CREATED             STATUS                      PORTS               NAMES
39548f30ce45        my_docker_conda     "/bin/bash -c 'bas..."    3 minutes ago       Exited (0) 3 minutes ago                             el
```

If we run `docker run` without any flags, your local terminal is attached to the
container. This enables you to see the output of `run_qc.sh`, but also disables
you from doing anything else in the meantime. We can start a container in
detached mode with the `-d` flag. Try this out and run `docker container ls`
to validate that the container is running.

By default, Docker keeps containers after they have exited. This can be
convenient for debugging or if you want to look at logs, but it also consumes
huge amounts of disk space. It's therefore a good idea to always run with
`--rm`, which will remove the container once it has exited.

If we want to enter a running container, there are two related commands we can
use, `docker attach` and `docker exec`. `docker attach` will attach local
standard input, output, and error streams to a running container. This can be
useful if your terminal closed down for some reason or if you started
a terminal in detached mode and changed your mind. `docker exec` can be used to
execute any command in a running container. It's typically used to peak in at
what is happening by opening up a new shell. Here we start the container in
detached mode and then start a new interactive shell so that we can see what
happens. If you use `ls` inside the container you can see how the script
generates file in the `data`, `intermediate` and `results` directories. Note
that you will be thrown out when the container exits, so you have to be quick.

```bash
docker run -d --rm --name my_container my_docker_conda
docker exec -it my_container /bin/bash
```

## Bind mounts

There are obviously some advantages to isolating and running your data analysis
in containers, but at some point you need to be able to interact with the host
system to actually deliver the results. This is done via bind mounts. When you
use a bind mount, a file or directory on the *host machine* is mounted into
a container. That way, when the container generates a file in such a directory
it will appear in the mounted directory on your host system.

!!! Tip
    Docker also has a more advanced way of data storage called
    [volumes](https://docs.docker.com/storage/volumes/). Volumes provide
    added flexibility and are independent of the host machine's filesystem
    having a specific directory structure available. They are particularly
    useful when you want to share data *between* containers.

Say that we are interested in getting the resulting html reports from FastQC in
our container. We can do this by mounting a directory called, say,
`fastqc_results` in your current directory to the `/course/results/fastqc`
directory in the container. Try this out by running:

```bash
docker run --rm -v $(pwd)/fastqc_results:/course/results/fastqc my_docker_conda
```

Here the `-v` flag to docker run specifies the bind mount in the form of
`directory/on/your/computer:/directory/inside/container`. `$(pwd)` simply
evaluates to the working directory on your computer.

Once the container finishes validate that it worked by opening one of the html
reports under `fastqc_results/`.

We can also use bind mounts for getting files into the container rather than
out. We've mainly been discussing Docker in the context of packaging an
analysis pipeline to allow someone else to reproduce its outcome. Another
application is as a kind of very powerful environment manager, similarly to how
we've used Conda before. If you've organized your work into projects, then you
can mount the whole project directory in a container and use the container as
the terminal for running stuff while still using your normal OS for editing
files and so on. Let's try this out by mounting our current directory and start
an interactive terminal. Note that this will override the `CMD` command, so we
won't start the analysis automatically when we start the container.

```bash
docker run -it --rm -v $(pwd):/course/ my_docker_conda /bin/bash
```

If you run `ls` you will see that all the files in the `docker` directory are
there. Now edit `run_qc.sh` **on your host system** to download, say, 12000
reads instead of 15000. Then rerun the analysis with `bash run_qc.sh`. Tada!
Validate that the resulting html reports look fine and then exit the container
with `exit`.

!!! Success "Quick recap"
    In this section we've learned:

    - How to use `docker run` for starting a container and how the flags `-d`
    and `--rm` work.
    - How to use `docker container ls` for displaying information about the
    containers.
    - How to use `docker attach` and `docker exec` to interact with running
    containers.
    - How to use bind mounts to share data between the container and the host system.
