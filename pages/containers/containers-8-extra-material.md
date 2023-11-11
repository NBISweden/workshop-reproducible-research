Containers can be large and complicated, but once you start using them regularly
you'll find that you start understand these complexities. There are lots of
different things you can do with images and containers in general, especially
when it comes to optimising build time or final image size. Here is some small
tips and tricks that you can be inspired from!

If you want to read more about containers in general you can check out these
resources:

* A "Get started with Docker" at the [Docker website](https://docs.docker.com/get-started/).
* An [early paper](https://arxiv.org/abs/1410.0846) on the subject of using
  Docker for reproducible research.

## Building for multiple platforms

With the newer ARM64 architectures introduced by Apple one often runs into the
problem of not having an architecture-native image to run with. This is
sometimes okay since the [Rosetta2](https://support.apple.com/en-us/HT211861)
software can emulate the old AMD64 architecture on newer ARM64 computers, but
results in a performance hit. One could just build for ARM64 using
`--platform=linux/arm64` instead, but then somebody who *doesn't* have the new
architecture can't run it. There is a way around this, however: *multi-platform
builds*. We can build for multiple platforms at the same time and push those to
*e.g.* DockerHub and anybody using those images will automatically pull the one
appropriate for their computer. Here's how to do it:

* Start by checking the available builders using `docker buildx ls`.

You should only see the default builder, which does not have access to
multi-platform builds. Let's create a new builder that *does* have access to it:

* Run the following: `docker buildx create --name mybuilder --driver
  docker-container --bootstrap`.

* Switch to using the new builder with `docker buildx use mybuilder` and check
  that it worked with `docker buildx ls`.

All that's needed now is to build and push the images! The following command
assumes that you have an account with `<username>` at DockerHub and you're
pushing the `<image>` image:

```bash
docker buildx build --platform linux/amd64,linux/arm64 -t <username>/<image>:latest --push .
```

* Execute the above command with your username and your image.

That's it! Now anybody who does *e.g.* `docker pull <username>/<image>` will get
an image appropriate for their architecture whether they are on AMD64 or ARM64!

> **An alias to `buildx`** <br>
> You can type `docker buildx install` to make the `docker build` into an alias
> for `docker buildx`, allowing you to run multi-platform builds using `docker
> build`. Use `docker buildx uninstall` to remove this alias.
