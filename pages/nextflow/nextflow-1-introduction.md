[Nextflow](https://www.nextflow.io/) is another *workflow management system*,
and is (along with Snakemake) one of the most common such systems within the
academic community. While Nextflow is based in the Groovy language, you don't
need to know how to code Groovy to be able to write good Nextflow workflows,
just like you don't need to know Python to write good Snakemake workflows.

Nextflow is built from the ground-up to be both portable, scalable, reproducible
and usable in a platform-agnostic sense. This means that any pipeline you write
in Nextflow can easily be run both locally on your laptop, a computer cluster or
a cloud service. You can also define the environment with which each task is
carried out, just like in Snakemake. Nextflow has a large community centered
around it, including the [nf-core](https://nf-co.re/) curated collection of high
quality pipelines used by *e.g.* the [National Genomics Infrastructure](https://ngisweden.scilifelab.se/).

#### Differences between Nextflow and Snakemake

There are some major differences in the philosophy behind Nextflow and
Snakemake. First of all, everything in Snakemake is a file, while Nextflow uses
data streams called *channels*. These channels allow you to not only work on
files, but also variables and other metadata. A common use of this is to work on
a file containing data from a sample along with its sample name, together, which
simplifies a lot of the coding. The flexibility of channels allow you to do many
different things in powerful ways, but a downside is that they can be
complicated when you start working with them.

Snakemake uses a "pull"-philosophy similar to its inspiring predecessor
[make](https://www.gnu.org/software/make/), meaning that you define a number of
rules with inputs and outputs and then ask for the specific result you want,
*i.e.* the final output files (usually defined in a `all` rule). Snakemake will
work backwards from the final outputs you desired and find whatever combination
of inputs and rules it needs to give them to you. This means that you always
know exactly which files are going to be created and manipulated in all steps of
the workflow even before it starts, which is a nice thing to know.

Nextflow works in the opposite way, *i.e.* with a "push"-philosophy: you define
a number of *processes* with inputs and outputs (equivalent to rules in
Snakemake), and then give the first inputs to Nextflow. It will run the first
process using those inputs, pass them to the second process, then the third, and
so on until you reach the final outputs of the workflow. This means that
Nextflow doesn't actually know exactly which files are going to be created and
manipualted during a run, which is both good and bad. The bad means that you
can't really do dry runs in Nextflow in the same simple manner as for Snakemake
(there is, however, something similar called [stubs](https://github.com/nextflow-io/nextflow/blob/master/docs/process.rst#stub)).
On the other hand, something very good about this is that Nextflow easily
handles variable inputs and outputs well, making dynamic analyses easy to work
with, *e.g.* processes where you don't know the exact number of output files. A
"push"-philosophy is also generally more intuitive for many people, although
this is mostly relevant when learning.

Another difference is that in Nextflow, every single process is run in its own,
hidden directory as an isolated environment; only the files relevant to that
process will be present in this directory. This greatly simplifies testing and
debugging, as you can always go into a process' directory and see exactly which
files it has access to and which code was executed. This general structure means
that you need to think less about full paths for all the workflow's files, as
the locations of all the files are fully taken care of by Nextflow - the only
thing you need to care about is the final output paths.

#### The aim of this tutorial

Will all that said, both Nextflow and Snakemake are excellent systems for
workflow management, and you can do basically everything in either: it is very
much up to personal preference and what you think is most important. We suggest
that you try both and get a feel for them, and then decide which you like the
most. The course material for Nextflow is, with this in mind, not as extensive
as that for Snakemake. We have recreated the MRSA workflow from Snakemake in
Nextflow and will, through it, give you an overview of how Nextflow does things.
The idea is that you should have a rough idea of both Snakemake and Nextflow, so
that you may continue in what manner you think suits you the best.
