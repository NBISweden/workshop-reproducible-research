<iframe id="iframepdf" src="../../../lectures/containers/containers.pdf" frameborder="0" width="640" height="480" allowfullscreen="true" mozallowfullscreen="true" webkitallowfullscreen="true"></iframe>

Container-based technologies are designed to make it easier to create, deploy,
and run applications by isolating them in self-contained software units (hence
their name). The idea is to package software and/or code together with
everything it needs (other packages it depends, various environment settings,
*etc.*) into one unit, *i.e.* a container. This way we can ensure that the
software or code functions in exactly the same way regardless of where it's
executed. Containers are in many ways similar to virtual machines but more
lightweight. Rather than starting up a whole new operating system, containers
can use the same kernel (usually Linux) as the system that they're running on.
This makes them much faster and smaller compared to virtual machines. While
this might sound a bit technical, actually using containers is quite smooth and
very powerful.

Containers have also proven to be a very good solution for packaging, running
and distributing scientific data analyses. Some applications of containers
relevant for reproducible research are:

* When publishing, package your analyses in a container image and let it
  accompany the article. This way interested readers can reproduce your analysis
  at the push of a button.
* Packaging your analysis in a container enables you to develop on *e.g.* your
  laptop and seamlessly move to cluster or cloud to run the actual analysis.
* Say that you are collaborating on a project and you are using Mac while your
  collaborator is using Windows. You can then set up a container image specific
  for your project to ensure that you are working in an identical environment.

One of the largest and most widely used container-based technologies is
*Docker*. Just as with Git, Docker was designed for software development but is
rapidly becoming widely used in scientific research. Another container-based
technology is *Singularity*, which was developed to work well in computer
cluster environments, such as Uppmax. We will cover both Docker and Singularity
in this course, but the focus will be be on the former (since that is the most
widely used and runs on all three operating systems).

This tutorial depends on files from the course GitHub repo. Take a look at the
[setup](pre-course-setup) for instructions on how to install Docker if you
haven't done so already, then open up a terminal and go to
`workshop-reproducible-research/tutorials/containers`.

!!! warning
    Docker images tend to take up quite a lot of space. In order to do all
    the exercises in this tutorial you need to have ~10 GB available.
