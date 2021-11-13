It is time to try to set up a project from scratch and use one (or several)
of the tools that we have covered during the course! This exercise is very
open-ended and you have free hands to try out whatever tool(s) you want. 

A fully reproducible project could be done like this:

1. Create a new Git repository for the project and connect it to GitHub;
   continuously commit changes to this repository as you work with your project
2. Add a `README.md` file, which should contain a description of the project
   and the required information on how to run all its analyses
3. Create a Conda `environment.yml` file with the required dependencies
4. Create a R Markdown or Jupyter notebook with your analyses
5. Create a Snakemake or Nextflow workflow to connect all the analyses
6. Create a Docker or Singularity image for your project

This is not a small task and may seem overwhelming! Don't worry if you feel
lost or if the task seems daunting. To get the most out of this exercise, 
think about which of the tools is most useful for your research project(s) 
right now, start working with that tool during this session and use the 
opportunity to ask questions. After the course, take one step at a time to
add more and more tools to your project and go back to the tutorials for 
help and inspiration.

> **Note** <br>
> Git, Conda and notebooks could be seen as the core tools to make a research 
> project reproducible, so we think that any of these tools is a good starting
> point. Whichever tool you choose to work with, we suggest to keep the analysis
> for this exercise short so that you have time to try out the tool while you
> have the opportunity to get help.

> **Reaching full reproducibility** <br>
> Getting to a fully reproducible project is not the easiest thing, but if
> you've done it once it'll be a lot easier! At that point, you'll find that not
> only do you now reach an extremely high scientific standard, but you'll also
> see how easy it is to work inside your projects! Adding additional analyses,
> changing parameters, fixing things for review, *etc.* will all become much
> easier for you. If you follow all of the points in the list at the top of this
> page, you'll get there in no time!


## Your own project

This is a great opportunity for you to try to implement the tools on one of
your current research projects. It is of course up to you which tool(s) to
include to make your research project reproducible.

> **Tip** <br>
> If your analysis project contains computationally intense steps it may be
> good to scale them down for the sake of the exercise. You might, for
> example, subset your raw data to only contain a minuscule part of its
> original size. You can then test your implementation on the subset and only
> run it on the whole dataset once everything works to your satisfaction.

## Alternative: student experience project

If you don't want to use a project you're currently working on we have
a suggestion for a small-scale project for you. The idea is to analyze
students' experiences at this Tools for Reproducible Research course. For this
you will use responses from students to the registration form from previous
course rounds. You can find the responses in the `workshop-reproducible-research/tutorials/data/`
directory. The goal is to:

1. Create a project workflow using Snakemake or Nextflow, containing a step to
   clean the files (making sure the rule/process is generalisable for different
   input data)
2. Create a plot of the student experiences in some interesting way,
   implemented in the same workflow as above or using R Markdown / Jupyter

