# Introduction

It is time to try to set up a project from scratch and use the different
tools that we have covered during the course together! This exercise is very
open-ended and you have free hands to try out a bit of what you want. But you
should aim to use what you've learned to do the following:

1. Create a new git repository for the project (either on BitBucket or GitHub)
2. Add a README file which should contain the required information on how to
   run the project
3. Create a Conda `environment.yml` file with the required dependencies
4. Create a R Markdown or Jupyter notebook to run your code
5. Alternatively, create a `Snakefile` to run your code as a workflow and use a `config.yml` file 
   to add settings to the workflow
6. Use git to continuously commit changes to the repository
7. Possibly make a Docker or Singularity image for your project

This is not a small task and may seem overwhelming! Don't worry if you feel
lost or if the task seems daunting. To get the most out of the exercise, take
one step at a time and go back to the previous tutorials for help and
inspiration. The goal is not necessarily for you to finish the whole exercise,
but to really think about each step and how it all fits together in practice.

> **Note** <br>
> We recommend to start with Git, Conda and a notebook, as they could be seen as
> the core tools to make a research project reproducible. We suggest to keep the
> analysis for this exercise short so that you have time to try out the
> different tools together while you have the opportunity to ask for help.

## Alternative 1: your own project

This is a great opportunity for you to try to implement these methods on one 
of your current research projects. It is of course up to you which tools to 
include in making your research project reproducible, but we suggest to aim 
for at least git and Conda. 

> **Tip** <br>
> If your analysis project contains computationally intense steps it may be
> good to scale them down for the sake of the exercise. You might, for
> example, subset your raw data to only contain a minuscule part of its
> original size. You can then test your implementation on the subset and only
> run it on the whole dataset once everything works to your satisfaction.

## Alternative 2: student experience project

If you don't want to use a project you're currently working on we have
a suggestion for a small-scale project for you. The idea is to analyze students'
experiences at this Reproducible Research course. For this you will use
responses from students to the registration form for previous course rounds. You
can find the responses in the `data/` directory in the project root. The goal
is to create a Snakemake workflow, which contains the following:

1. A rule that cleans the files (making use of `wildcards` so that the same rule
   can be run on each file)

2. Plot the student experiences in some interesting way

> **Attention** <br>
> Remember to: <br>
>
> * Keep everything versioned controlled with `git` <br>
> * Add information to the `README` file so others know how to re-run 
>   the project <br>
> * Add required software to the Conda `environment.yml` file <br>
