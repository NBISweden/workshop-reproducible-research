Docker is a tool designed to make it easier to create, deploy, and run
applications by isolating them in "containers". The idea is to package your
code together with everything it needs (other packages it depends on, various
environment settings, *etc.*) into one unit, *i.e.* a container. This way we
can ensure that the code results in exactly the same output regardless of where
it's executed. Containers are in many ways similar to virtual machines but more
lightweight. Rather than starting up a whole new OS, Docker containers can use
the same Linux kernel as the system that they're running on. This makes them
much faster and smaller compared to virtual machines. While this might sound
a bit technical, actually using Docker is quite easy, fun and very powerful.

Just as with Git, Docker was designed for software development but is rapidly
becoming widely used in scientific research. Say that you are building a web
application: you could then run the web server in one container and the
database in another, thereby reducing the risk of one system affecting the
other in unpredictable ways. Docker containers have also proven to be a very
good solution for packaging, running and distributing scientific data analyses.
Some applications relevant for reproducible research can be:

* When publishing, package your whole analysis pipeline in a Docker image and
  let it accompany the article. This way interested readers can reproduce your
  analysis at the push of a button.
* Packaging your analysis in a Docker container enables you to develop on
  *e.g.* your laptop and then seamlessly move to cluster or cloud to run the
  actual analysis.
* Say that you are collaborating on a project and you are using Mac while your
  collaborator is using Windows. You can then set up a Docker image specific
  for your project to ensure that you are working in an identical environment.

All of this might sound a bit abstract so far, but it'll become more clear
during the exercises. If you want to read more you can check out these
resources:

* A "Get started with Docker" at [docker.com](https://docs.docker.com/get-started/).
* An [early paper](https://arxiv.org/abs/1410.0846) on the subject of using
  Docker for reproducible research.

This tutorial depends on files from the course GitHub repo. Take a look at the
[setup](pre-course-setup) for instructions on how to install Docker if you 
haven't done so already, then open up a terminal and go to 
`workshop-reproducible-research/docker`.

> **Attention!** <br>
> Docker images tend to take up quite a lot of space. In order to do all
> the exercises in this tutorial you need to have ~10 GB available.
