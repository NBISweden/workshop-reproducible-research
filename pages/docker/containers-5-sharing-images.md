There would be little point in going through all the trouble of making your
analyses reproducible if you can't distribute them to others. Luckily, sharing
Docker containers is extremely easy. The most common way is to use
[Dockerhub](https://hub.docker.com). Dockerhub lets you host unlimited public
images and one private image for free, after that they charge a small fee. If
you want to try it out here is how to do it:

1. Register for an account on [Dockerhub](https://hub.docker.com).

2. Use `docker login -u your_dockerhub_id` to login to the Dockerhub registry.

3. When you build an image, tag it with `-t your_dockerhub_id/image_name`,
   rather than just `image_name`.

4. Once the image has been built, upload it to Dockerhub with `docker push
   your_dockerhub_id/image_name`.

5. If another user runs `docker run your_dockerhub_id/image_name` the image
   will automatically be retrieved from Dockerhub. You can use `docker pull`
   for downloading without running.

If you want to refer to a Docker image in for example a publication, it's very
important that it's the correct version of the image. You can do this by adding
a tag to the name like this `docker build -t
your_dockerhub_id/image_name:tag_name`.

> **Tip** <br>
> On Dockerhub it is also possible to link to your Bitbucket or GitHub
> account and select repositories from which you want to automatically build
> and distribute Docker images. The Dockerhub servers will then build an
> image from the Dockerfile in your repository and make it available for
> download using `docker pull`. That way, you don't have to bother manually
> building and pushing using `docker push`. The GitHub repository for this
> course is linked to Dockerhub and the Docker images are built automatically
> from `Dockerfile` and `Dockerfile_slim`, triggered by changes made to the
> GitHub repository. You can take a look at the course on Dockerhub
> [here](https://hub.docker.com/r/nbisweden/workshop-reproducible-research).

> **Quick recap** <br>
> In this section we've learned:
>
> - How to use [DockerHub](https://hub.docker.com) to share Docker images
