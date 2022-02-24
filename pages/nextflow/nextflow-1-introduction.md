[Nextflow](https://www.nextflow.io/) is another *workflow management system*,
and is (along with Snakemake) one of the most common such systems within the
academic community.

Nextflow is built from the ground-up to be portable, scalable, reproducible and
usable in a platform-agnostic sense. This means that any pipeline you write in
Nextflow can be configured to run locally on your laptop, a computer cluster or
a cloud service (as long as your architecture has the necessary compute
resources). You can also define the compute environment in which each task is
carried out, just like in Snakemake. Nextflow has a large community centered
around it, including the [nf-core](https://nf-co.re/) curated collection of
high quality pipelines used by *e.g.* the [National Genomics Infrastructure](https://ngisweden.scilifelab.se/).

# Differences between Nextflow and Snakemake

There are several major differences between Snakemake and Nextflow, both in how
they work and in their underlying philosophies, summarised by the following
table:

<table class="table table-hover table-condensed" border=1; style="width:600px; margin-left:auto; margin-right:auto;">
    <thead style="background-color:#DAE7F1">
        <tr>
            <td style="padding:5px"> <font size="3"></td>
            <td style="padding:5px"> <font size="3"><b> Snakemake </b> </td>
            <td style="padding:5px"> <font size="3"><b> Nextflow </b> </td>
        </tr>
    </thead>
    <tr>
        <td style="padding:5px"> <font size="3"> <b>Language</b> </td>
        <td style="padding:5px"> <font size="3"> Python </td>
        <td style="padding:5px"> <font size="3"> Groovy </td>
    </tr>
    <tr>
        <td style="padding:5px"> <font size="3"> <b>Data</b> </td>
        <td style="padding:5px"> <font size="3"> Everything is a file </td>
        <td style="padding:5px"> <font size="3"> Can use both files and values </td>
    </tr>
    <tr>
        <td style="padding:5px"> <font size="3"> <b>Execution</b> </td>
        <td style="padding:5px"> <font size="3"> Working directory </td>
        <td style="padding:5px"> <font size="3"> Each task is isolated in it's own directory </td>
    </tr>
    <tr>
        <td style="padding:5px"> <font size="3"> <b>Philosophy</b> </td>
        <td style="padding:5px"> <font size="3"> "Pull" </td>
        <td style="padding:5px"> <font size="3"> "Push" </td>
    </tr>
    <tr>
        <td style="padding:5px"> <font size="3"> <b>Dry runs</b>  </td>
        <td style="padding:5px"> <font size="3"> Yes </td>
        <td style="padding:5px"> <font size="3"> No </td>
    </tr>
    <tr>
        <td style="padding:5px"> <font size="3"> <b>Track code changes</b> </td>
        <td style="padding:5px"> <font size="3"> No </td>
        <td style="padding:5px"> <font size="3"> Yes </td>
    </tr>
</table>

Starting from the top with the most obvious and perhaps superficial difference:
**language**. Snakemake is, as you know, based on Python, whereas Nextflow is
based on Groovy (which is a superset of the Java language). You don't need to
know Groovy to be able to use Nextflow though, just as you don't really need to
know Python to use Snakemake.

Moving on to **data**: Input and output in Snakemake are files, whereas
Nextflow uses objects, which can be files, environment variables, or arbitrary
data structures, transmitted through so-called *channels*. Channels are
asynchronous first-in-first-out streams of data that connect a workflow's
various inputs and outputs. A common use-case is, for example, to define
a channel which passes objects containing both sample data files and their
corresponding sample names, which can simplify coding. Nextflow also defines
channel operators; functions that allow you to manipulate channel contents in
powerful ways, although they can be tricky to use when you first start working
with them.

In Snakemake, the entire workflow and each rule is **executed** in the working
directory, while Nextflow executes each individual *task* (an instance of
a *process* - the equivalent of a rule in Snakemake) within an isolated
environment in a directory of its own. This greatly simplifies testing and
debugging, as you can always go into a process' directory and see exactly which
files it has access to and which code was executed. This general structure
means that you need to think less about full paths for all the workflow's in-
and output files, as the locations of all the files are fully taken care of by
Nextflow - the only thing you need to care about are the final output file
paths.

Snakemake uses a "pull"-**philosophy** similar to its inspiring predecessor
[make](https://www.gnu.org/software/make/), meaning that you define a number of
rules with inputs and outputs and then ask for the specific result you want,
*i.e.* the final output files (usually defined in an `all` rule). Snakemake will
work backwards from the final outputs you desired and find whatever combination
of inputs and rules it needs to give them to you. This means that you always
know exactly which files are going to be created and manipulated in all steps of
the workflow even before it is executed, which is a nice thing to know. Nextflow
works in the opposite way, *i.e.* with a "push"-philosophy: you define a number
of processes with inputs and outputs, and then give the first inputs to
Nextflow. It will run the first process using those inputs, pass them to the
second process, then the third, and so on until it reaches the final outputs of
the workflow. This means you don't define file paths to each process'
input/output definitions like you do in Snakemake, only which files you want in
the end. This can potentially remove some of the pitfalls and issues sometimes
seen with *e.g.* wildcards in Snakemake, but it also front-loads some of the
complexity to channel creation instead.

The philosophy above means that Nextflow doesn't know exactly which files
are going to be created and manipulated during a run, which is both good and
bad. The bad means that you can't really do **dry runs** in Nextflow in the same
simple manner as in Snakemake (there is, however, something similar called
[stubs](https://github.com/nextflow-io/nextflow/blob/master/docs/process.rst#stub)).
On the other hand, something good about this is that Nextflow handles variable
inputs and outputs very well, making dynamic analyses easy to work with, *e.g.*
processes where you don't know the exact number of output files.

Lastly, both Snakemake and Nextflow can automatically determine which rules or
processes need to be re-run when something has changed, but they do it in
slightly different ways. Snakemake only checks if any of the input files are
newer than the output files, while Nextflow also **tracks code updates, changes
of software environment, and changes in input values**. This means that if you
update a script that is run on some unchanged data in Nextflow, it will re-run
the corresponding process automatically; the same is not true for Snakemake,
where you need to specify that you want to re-run the workflow from that
specific rule (*i.e.* using `-R <rule>`).

# The aim of this tutorial

With all that said, both Nextflow and Snakemake are excellent systems for
workflow management, and you can do basically everything in either: your choice
is very much up to your personal preference and what you think is most
important. We suggest that you try both and get a feel for them, and then
decide which you like the most. The idea is that you should have a rough idea
of both Snakemake and Nextflow after the course, so that you may continue in
what manner you think suits you the best.

This tutorial depends on files from the course GitHub repo. Take a look at the
[setup](pre-course-setup) for instructions on how to set it up if you haven't
done so already, then open up a terminal and go to `workshop-reproducible-research/tutorials/nextflow`
and activate your `nextflow-env` Conda environment.
