# Making a project reproducible

It is time to try to setup a project from scratch and use the different 
tools that we have covered during the course together! 

This exercise if very open-ended and you have free hands to try out a 
bit what you want. 

## Option 1
One option is to try to implement these methods on 
one of your current projects. It is up to you what tools to include in 
making your project reproducible, but aim for at least including git and conda.

!!! tip
    If your analysis project contains 
    computationally intense steps it may be good to scale down them for the sake of the exercise. 


## Option 2
If you don't want to use a project you're currently working on we have 
a suggestion for a small-scale project for you.

The idea is to analyze student experience for this Reproducible Research
course. For this you will use responses from students to the registration 
form for the course. Below you'll find links to **csv** format files
with answers from 3 course instances:

* 2018-11:https://docs.google.com/spreadsheets/d/1yLcJL-rIAO51wWCPrAdSqZvCJswTqTSt4cFFe_eTjlQ/export?format=csv
* 2019-05:https://docs.google.com/spreadsheets/d/1mBp857raqQk32xGnQHd6Ys8oZALgf6KaFehfdwqM53s/export?format=csv
* 2019-11:https://docs.google.com/spreadsheets/d/1aLGpS9WKvmYRnsdmvvgX_4j9hyjzJdJCkkQdqWq-uvw/export?format=csv

## Objectives

1. Create a new git repository for the project (either on BitBucket or GitHub)
2. Initialize a README file which should contain the required information on how to run the project
3. Create a conda `environment.yml` file with the required dependencies
4. Create a snakemake file 

??? note Click to show a script for cleaning column names