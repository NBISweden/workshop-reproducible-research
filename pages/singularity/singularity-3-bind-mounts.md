# Bind mounts

In the previous section we saw how Singularity differs from Docker in terms of
images being stored in stand-alone files and much of the host filesystem being
mounted in the container. We will now explore this further.

Similarly to the `-v` flag for Docker we can use `-B` or `--bind` to bind mount
directories into the container. For example, to mount a directory into the
`/mnt/` directory in the container one would do:

```bash
singularity shell -B /path/to/dir:/mnt ubuntu_latest.sif
```

You can try this for example by mounting the `conda/` tutorial directory to
`/mnt`:

```bash
singularity shell -B ../conda:/mnt ubuntu_latest.sif
```

In the container, to see that the bind mounting worked, run *e.g.*:

```bash
ls /mnt/code
```

Now, this was not really necessary since `conda/` would have been available to
us anyway since it most likely is a sub-directory under your `$HOME`, but it
illustrates the capabilities to get files from the host system into the
container when needed. Note also that if you run Singularity on say an HPC
cluster, the system admins may have enabled additional default directories that
are bind mounted automatically.

> **Quick recap** <br>
> In this section we covered how to bind mount specific directories using `-B`.
