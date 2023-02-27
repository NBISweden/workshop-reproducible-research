There would be little point in going through all the trouble of making your
analyses reproducible if you can't distribute them to others. Luckily, sharing
Docker containers is extremely easy, and can be done in several ways. One of the
more common ways to share Docker images is through container *registries* and
*repositories*.

For example, a Docker registry is a service that stores Docker images, which
could be hosted by a third party, publicly or privately. One of the most common
registries is [Docker Hub](https://docs.docker.com/docker-hub/), which is a
registry hosted by Docker itself. A repository, on the other hand, is a
collection of container images with the same name but different tags, *i.e.*
versions. For example, `ubuntu:latest` or `ubuntu:20.04`. Repositories are
stored in registries.

!!! Note
    Remember that we now have some clashing nomenclature between Git repositories
    (which we covered in the Git tutorial) and container repositories, so be aware
    of which one you're talking about!

There are many registries out there, but here are some that might be of interest
to you who are taking this course:

* [Docker Hub](https://docs.docker.com/docker-hub/)
* [Quay](https://quay.io/)
* [Biocontainers](https://biocontainers.pro/#/registry)
* [Rocker](https://www.rocker-project.org/images/)
* [Jupyter containers](https://jupyter-docker-stacks.readthedocs.io/en/latest)

The most common registry is probably Docker Hub, which lets you host unlimited
public images and one private image for free (after which they charge a small
fee). Let's see how it's done!

1. Register for an account on [Docker Hub](https://hub.docker.com).

2. Use `docker login -u your_dockerhub_id` to login to the Docker Hub registry.

3. When you build an image, tag it with `-t your_dockerhub_id/image_name`,
   rather than just `image_name`.

4. Once the image has been built, upload it to Docker Hub with `docker push
   your_dockerhub_id/image_name`.

5. If another user runs `docker run your_dockerhub_id/image_name` the image
   will automatically be retrieved from Docker Hub. You can use `docker pull`
   for downloading without running.

If you want to refer to a Docker image in for example a publication, it's very
important that it's the correct version of the image. You can do this by adding
a tag to the name like this `docker build -t
your_dockerhub_id/image_name:tag_name`.

!!! Tip
    On Docker Hub it is also possible to link to your Bitbucket or GitHub
    account and select repositories from which you want to automatically build
    and distribute Docker images. The Docker Hub servers will then build an
    image from the Dockerfile in your Git repository and make it available for
    download using `docker pull`. That way, you don't have to bother manually
    building and pushing using `docker push`. The GitHub repository for this
    course is linked to Docker Hub and the Docker images are built automatically
    from `Dockerfile` and `Dockerfile_slim`, triggered by changes made to the
    GitHub repository. You can take a look at the course on Docker Hub
    [here](https://hub.docker.com/r/nbisweden/workshop-reproducible-research).

!!! Success "Quick recap"
    In this section we've learned:

    - How container registries and repositories work
    - How to use Docker Hub to share Docker images
