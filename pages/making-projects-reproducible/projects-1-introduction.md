It is time to try to set up a project from scratch and use one (or several)
of the tools that we have covered during the course! This exercise is very
open-ended and you have free hands to try out whatever tool(s) you want. 

A fully reproducible project could look like that:

1. Create a new git repository for the project (either on BitBucket or GitHub)
2. Add a README file, which should contain the required information on how to
   run the project
3. Create a Conda `environment.yml` file with the required dependencies
4. Create a R Markdown or Jupyter notebook to run your code
5. Alternatively, create a Snakemake or Nextflow file to run your code as 
   a workflow
6. Use git to continuously commit changes to the repository
7. Possibly make a Docker or Singularity image for your project

This is not a small task and may seem overwhelming! Don't worry if you feel
lost or if the task seems daunting. To get the most out of this exercise, 
think about which of the tools is most useful for your research project(s) 
right now, start working with that tool during this session and use the 
opportunity to ask questions. After the course, take one step at a time to
add more and more tools to your project and go back to the tutorials for 
help and inspiration.

> **Note** <br>
> Git, Conda and notebooks could be seen as the core tools to make a research 
> project reproducible, so we think that any of these tools is a good starting point. 
> Whichever tool you choose to work with, we suggest to keep the analysis for 
> this exercise short so that you have time to try out the tool while you have 
> the opportunity to get help.

## Your own project

This is a great opportunity for you to try to implement one of the tools on one 
of your current research projects. It is of course up to you which tool(s) to 
include to make your research project reproducible.

> **Tip** <br>
> If your analysis project contains computationally intense steps it may be
> good to scale them down for the sake of the exercise. You might, for
> example, subset your raw data to only contain a minuscule part of its
> original size. You can then test your implementation on the subset and only
> run it on the whole dataset once everything works to your satisfaction.

## Alternative: student experience project

If you don't want to use a project you're currently working on we have
a suggestion for a small-scale project for you. The idea is to analyze students'
experiences at this Reproducible Research course. For this you will use
responses from students to the registration form for previous course rounds. You
can find the responses in the `workshop-reproducible-research/tutorials/data/`
directory. 

The goal is to create 

1. a Snakemake workflow, which contains a rule that cleans the files (making 
   use of `wildcards` so that the same rule can be run on each file). Alternatively,
   Nextflow can be used to run the file cleaning.

2. a plot of the student experiences in some interesting way, implemented in
   the same pipeline as above or using a notebook.


> **Note** <br>
> To reach full reproducibility:
>
> * Keep everything versioned controlled with `git`
> * Add information to the `README` file so others know how to re-run 
>   the project
> * Add required software to the Conda `environment.yml` file
