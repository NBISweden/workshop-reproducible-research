# Introduction to Snakemake

A workflow management system (WMS) is a piece of software that sets up, performs and monitors a defined sequence of computational tasks (i.e. "a workflow"). Snakemake is a WMS that was developed in the bioinformatics community, and as such it has some features that make it particularly well suited for creating reproducible and scalable data analyses.

* The language you use to formulate your workflows is based on Python; a language with a strong standing in academia. Users are not required to know how to code in Python to work efficently with Snakemake, they can still take advantage of Python's clean syntax and large number of libraries.

* Workflows can easily be scaled from your desktop to server, cluster, grid or cloud environments. This makes it possible to develop a workflow on your laptop, maybe using only a small subset of your data, and then run the real analysis on a cluster.

* Snakemake has several features for defining the environment which each task is carried out in. This is important in bioinformatics, where workflows often involve running a large number of small third-party tools.

* Snakemake is primarily intended to work on _files_ (rather than for example streams, reading/writing from databases or passing variables in memory). This fits well with many fields of bioinformatics, notably next-generation sequencing, which often involve computationally expensive operations on large files. It's also a good fit for a research setting and exploratory analyses, where the exact specifications of the final workflow aren't always known in the beginning of a project.

* Lastly, a WFM is a very important tool for making your analyses reproducible. By keeping track of when each file was generated, and by which operation, it is possible to ensure that there is a consistent "paper trail" from raw data to final results. Snakemake also has features which allow you to package and distribute the workflow, and any files it involves, once it's done.

## Tell me more
* The Snakemake documentation is available on [readthedocs](https://snakemake.readthedocs.io/en/stable/#).
* Here is a great [tutorial](https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html#tutorial).
* If you have questions, check out [stack overflow](https://stackoverflow.com/questions/tagged/snakemake).
* To discuss with other users, there is a [Google group](https://groups.google.com/forum/#!forum/snakemake).

# Set up
This tutorial depends on files from the course BitBucket repo. Take a look at the [intro](index) for instructions on how to set it up if you haven't done so already. Then open up a terminal and go to `reproducible_research_course/snakemake`.

## Set up the environment for the workflow
If you have done the [Conda tutorial](conda) you should know how to define an environment and install packages using Conda and an `environment.yml` file. Here we will use Snakemake as well as some other programs, all of which are available in the Bioconda channel. If you look in the current directory you will see an `environment.yml` file which specifies an environment containing FastQC, Bowtie2 and SRA Tools (identical to the output from the Conda tutorial). Add the following programs to the file and save it.

```yaml
# The workflow manager
- snakemake=4.3.0
# For aggregating output from FastQC
- multiqc=1.3
# For mapping reads to a genome
- bowtie2=2.3
# For sorting the output from Bowtie2
- samtools=1.6
# For generating a count table for further analysis
- htseq=0.9
# For visualizing workflows
- graphviz=2.38.0
- xorg-libxrender
- xorg-libxpm
```

Now install the new environment and activate it.

```bash
conda env create -n snakemake_exercise -f environment.yml
source activate snakemake_exercise
```

Check that Snakemake is installed correctly, for example by executing `bash snakemake --help`. This should output a list of available Snakemake settings. If you get `bash: snakemake: command not found` then you need to go back and ensure that the Conda steps were successful.

(If you don't want to use Conda for some reason you can also install Snakemake with `pip3 install snakemake`.)

# Practical exercise
## Toy example
In this part of the tutorial we will create a very simple workflow from scratch in order to show the fundamentals of how Snakemake works. The workflow will take two files as inputs, `a.txt` and `b.txt`, and the purpose is to convert the text in the files to upper case and then to concatenate them.

First make an empty file named `Snakefile` which will contain the workflow, and the two input files containing some arbitary text.

```bash
touch Snakefile
echo "This is a.txt" > a.txt
echo "This is b.txt" > b.txt
```

Then open `Snakefile` in your favorite text editor. A Snakemake workflow is based on rules which take some file(s) as input, performs some type of operation on them, and generate some file(s) as outputs. Here is a very simple rule that takes `a.txt` as an input and produces `a.upper.txt` as an output. Copy this rule to your `Snakefile` and save it.

```python
rule convert_to_upper_case:
    input:
        "a.txt"
    output:
        "a.upper.txt"
    shell:
        """
        tr [a-z] [A-Z] < {input} > {output}
        """
```

A rule has a name, here it's `convert_to_upper_case`. Make an effort to name your rules in way that makes it easy to understand the purpose of the rule, as rule names are one of the main ways to interact with the workflow. Intendation is important in the Snakemake language so make sure that you have tabs or spaces before `input`/`output`/`shell` and their respective subsections.

The `shell` section contains the shell commands that will convert the text in the input file to upper case and send it to the output file. In the shell command string, we can refer to elements of the rule via curly brackets. Here, we refer to the output file by specifying `{output}` and to the input file by specifying `{input}`. If you're not very familiar with Bash, this particular command can be read like "send the contents of `a.txt` to the program `tr`, which will convert all characters in the set [a-z] to the corresponding character in the set [A-Z], and then send the output to `a.upper.txt`".

Now let's run our first Snakemake workflow. When a workflow is executed Snakemake tries to generate given target files. Target files can be specified via the command line (or, as you will see later, in several other ways). Here we ask Snakemake to make the file `a.upper.txt`. It's good practice to first run with the flag `-n`(or `--dry-run`), which will show what Snakemake plans to do without actually running anything. You can also use the flag `-p` for showing the shell commands that it will execute and the flag `-r` for showing the reason for running a specific rule. `snakemake --help` will show you all available flags.

```no-highlight
$ snakemake -n -r -p a.upper.txt

rule convert_to_upper_case:
    input: a.txt
    output: a.upper.txt
    jobid: 0
    reason: Missing output files: a.upper.txt


        tr [a-z] [A-Z] < a.txt > a.upper.txt

Job counts:
        count   jobs
        1       convert_to_upper_case
        1
```

You can see that Snakemake plans to run 1 job: the rule `convert_to_upper_case` with `a.txt` as input and `a.upper.txt` as output. The reason for doing this is that it's missing the file `a.upper.txt`. Now execute the workflow without the `-n` flag and check that the content of `a.upper.txt` is as expected. Then try running the same command again. What do you see? It turns out that Snakemake only re-runs jobs if **one of the input files is newer than one of the output files, or if one of the input files will be updated by another job**. This is how Snakemake ensures that everything in the workflow is up to date. We will get back to this shortly.

What if we ask Snakemake to generate the file b.upper.txt?

```no-highlight
$ snakemake -n -r -p b.upper.txt
MissingRuleException:
No rule to produce b.upper.txt (if you use input functions make sure that they do not raise unexpected exceptions)
```

That didn't work well. We could copy the rule to make a similar one for `b.txt`, but that would be a bit cumbersome. Here is where named wildcards come in; one of the most powerful features of Snakemake. Simply change the input from `input: "a.txt"` to `input: "{some_name}.txt"` and the output to `output: "{some_name}.upper.txt"`. Now try asking for `b.upper.txt` again.

Tada! What happens here is that Snakemake looks at all the rules it has available (actually only one in this case) and tries to assign values to all wildcards so that the targeted files can be generated. In this was it was quite simple, you can see that it says that `wildcards: some_name=b`, but for large workflows and multiple wildcards it can get much more complex. Named wildcards is what enables a workflow (or single rules) to be efficiently generalized and reused between projects or shared between people.

It seems we have the first part of our workflow working, now it's time to make the second rule for concatenating the outputs from `convert_to_upper_case`. The rule structure will be similar, the only difference is that here we have two inputs instead of one. This can be expressed in two ways, either with named inputs like this:

```python
input:
    firstFile="...",
    secondFile="..."
shell:
    """
    some_function {input.firstFile} {input.secondFile}
    """
```

Or with indexes like this:

```python
input:
    "...",
    "..."
shell:
    """
    some_function {input[0]} {input[1]}
    """
```

Now try to construct this rule yourself and name it `concatenate_files`. The syntax for concatenating two files in Bash is `cat first_file second_file > output_file`. Call the output `c.txt`. Run the workflow in Snakemake and validate that the output looks as expected.

Wouldn't it be nice if our workflow could be used for _any_ files, not just `a.txt` and `b.txt`? We can achieve this by using named wildcards (or in other ways as we will discuss later). As we've discussed, Snakemake looks at all the rules it has available and tries to assign values to all wildcards so that the targeted files can be generated. We therefore have to name the output file in a way so that it also contains information about which input files it should be based on. Try to figure out how to do this yourself. If you're stuck you can look at the spoiler below, but spend some time on it before you look.

??? note "Click to see how to implement `concatenate_files`"
    ```python
    rule concatenate_files:
        input:
            "{first}.upper.txt",
            "{second}.upper.txt"
        output:
            "{first}_{second}.txt"
        shell:
            """
            cat {input[0]} {input[1]} > {output}
            """
    ```

We can now control which input files to use by the name of the file we ask Snakemake to generate.

```no-highlight
$ snakemake a_b.txt
Provided cores: 1
Rules claiming more threads will be scaled down.
Job counts:
        count   jobs
        1       concatenate_files
        2       convert_to_upper_case
        3

rule convert_to_upper_case:
    input: b.txt
    output: b.upper.txt
    jobid: 1
    wildcards: some_name=b

Finished job 1.

1 of 3 steps (33%) done

rule convert_to_upper_case:
    input: a.txt
    output: a.upper.txt
    jobid: 2
    wildcards: some_name=a

Finished job 2.
2 of 3 steps (67%) done

rule concatenate_files:
    input: a.upper.txt, b.upper.txt
    output: a_b.txt
    jobid: 0
    wildcards: second=b, first=a

Finished job 0.
3 of 3 steps (100%) done
```

Neat! All that we've done so far could quite easily be done in a simple shell script that takes the input files as parameters. Let's now take a look at some of the features where a WFM like Snakemake really adds value compared to a more straight forward approach. One such feature is the possibility to visualize your workflow. Snakemake can generate two types of graphs, one that show how the rules are connected and one that shows how the jobs (i.e. an execution of a rule with some given inputs/outputs/settings) are connected. First we look at the rule graph. The following command will generate a rule graph in the dot language and pipe it to the program `dot`, which in turn will save a visualisation of the graph as a png file (if you're having troubles displaying png files you could use svg or jpg instead).

```bash
snakemake --rulegraph a_b.txt | dot -Tpng > rulegraph.png
```

![alt text](rulegraph.svg)

That looks simple enough, the output from the rule `convert_to_upper_case` will be used as input to the rule `concatenate_files`. For a more typical bioinformatics project it can look something like this when you include all the rules from processing of the raw data to generating figures for the paper.

![alt text](rulegraph_complex.svg)

While saying that it's easy to read might be a bit of a stretch, it definitely gives you a better overview of the project than you would have without a WFM.

The second type of graph is based on the jobs, and looks like this for our little workflow (use `--dag` instead of `--rulegraph`).

```bash
snakemake --dag a_b.txt | dot -Tpng > jobgraph.png
```

![alt text](jobgraph.svg)

The main difference here is that now each node is a job instead of a rule. You can see that the wildcards used in each job are also displayed. Another difference is the dotted lines around the nodes. A dotted line is Snakemake's way of indicating that this rule doesn't need to be re-run in order to generate `a_b.txt`. Validate this by running `snakemake -n -r a_b.txt` and it should say that there is nothing to be done.

We've discussed before that one of the main purposes of using a WFM is that it automatically makes sure that everything is up to date. This is done by recursively checking that outputs are always newer than inputs for all the rules involved in the generation of you target files. Now try to change the content of `a.txt` to some other text and save it. What do you think will happen if you run `snakemake -n -r a_b.txt` again?

??? note "Click to see output"
    ```no-highlight
    $ snakemake -n -r a_b.txt

    rule convert_to_upper_case:
        input: a.txt
        output: a.upper.txt
        jobid: 2
        reason: Updated input files: a.txt
        wildcards: some_name=a



    rule concatenate_files:
        input: a.upper.txt, b.upper.txt
        output: a_b.txt
        jobid: 0
        reason: Input files updated by another job: a.upper.txt
        wildcards: first=a, second=b

    Job counts:
            count   jobs
            1       concatenate_files
            1       convert_to_upper_case
            2
    ```

Were you correct? Also generate the job graph and compare to the one generated above. What's the difference? Now rerun without `-n` and validate that `a_b.txt` contains the new text. Note that Snakemake doesn't look at the contents of files when trying to determine what has changed, only at the timestamp for when they were last modified.

We've seen that Snakemake keeps track of if files in the workflow have changed, and automatically makes sure that any results depending on such files are regenerated. What about if the rules themselves are changed? It turns out that there are multiple ways to do this, but the most straight forward is to manually specify that you want to rerun a rule (and thereby also all the steps between that rule and your target). Let's say that we want to modify the rule `concatenate_files` to also include which files were concatenated.

```python
rule concatenate_files:
    input:
        "{first}.upper.txt",
        "{second}.upper.txt"
    output:
        "{first}_{second}.txt"
    shell:
        """
        echo 'Concatenating {input}' | cat - {input[0]} {input[1]} > {output}
        """
```

It's not really important for the exercise, but the shell command used here first outputs "Concatenating " followed by a space delimited list of the files in `input`. This string is then sent to the program `cat` where it's concatenated with `input[0]` and `input[1]` (the parameter `-` means that it should read from standard input). Lastly, the output from `cat` is sent to `{output}`.

If you now run the workflow as before you should get "Nothing to be done" because no files involved in the workflow have been changed. Instead we have to force Snakemake to rerun the rule by using the `-R`flag. Let's try a dry-run.

```bash
snakemake a_b.txt -r -n -R concatenate_files
```

Note that the reason for the job is now "Forced execution". You can target files as well as rules, so you would get the same result with `-R a_b.txt`. Whenever you've made changes to a rule that will affect the output it's good practice to force re-execution like this. Still, there can be situations where you don't know if any rules have been changed. Maybe several people collaborate on the same workflow but are using it on different files for example. Snakemake keeps track of how all files were generated (when, by which rule, which version of the rule, and by which commands). You can export this information to a tab-delimited file like this:

```bash
snakemake a_b.txt -D > summary.tsv
```

The contents of `summary.tsv`is shown in the table below. You may want to open it in a spreadsheet like Microsoft Excel or Numbers to make it easier to read.


| output_file   | date   | rule   | version   | log-file(s)   |  input-file(s)   | shellcmd   | status   | plan   |
| --------------| ------ | ------ | --------- | ------------- | ---------------- | ---------- | -------- | ------ |
| a_b.txt | Thu Nov 16 12:03:11 2017 | concatenate_files | - |  | a.upper.txt,b.upper.txt | cat a.upper.txt b.upper.txt > a_b.txt | rule implementation changed | no update |
| a.upper.txt | Thu Nov 16 12:03:11 2017 | convert_to_upper_case | - |  | a.txt | tr [a-z] [A-Z] < a.txt > a.upper.txt | ok | no update |
| b.upper.txt | Thu Nov 16 12:03:11 2017 | convert_to_upper_case | - |  | b.txt | tr [a-z] [A-Z] < b.txt > b.upper.txt | ok | no update |


You can see in the second last column that the rule implementation for a_b.txt has changed. The last column shows if Snakemake plans to regenerate the files when it's next executed. None of the files will be regenerated because Snakemake doesn't regenerate files by default if the rule implementation changes. From a reproducibilty viewpoint maybe it would be better if this was done automatically, but it would be very computationally expensive and cumbersome if you had to rerun your whole workflow every time you fix a spelling mistake in a comment somewhere. So, it's up to us to look at the summary table and rerun things as needed. You can get a list of the files for which the rule implementation has changed, and then force Snakemake to regenerate these files with the `-R` flag.

```bash
snakemake a_b.txt -R `snakemake a_b.txt --list-code-changes`
```

Clever, right? There are a bunch of these flags that you can see with `--help`. Run with the `-D` flag again to make sure that the summary table now looks like expected.

You might wonder where Snakemake keeps track of all these things? It stores all information in a hidden subdirectory called `.snakemake`. This is convenient since it's easy to delete if you don't need it anymore and everything is contained in the project directory. Just be sure to add it to `.gitignore` so that you don't end up tracking it.

By now you should be familiar with the basic functionality of Snakemake, and you can build advanced workflows with only the features we have discussed here. There's a lot we haven't covered though, in particular when it comes to making your workflow more flexible and reusable. In the following section we will start with a workflow that is functional but not very flexible. We will then gradually improve on it and at the same time showcase some Snakemake features we haven't discussed yet. Note that this can get a little complex at times, so if you felt that this section was a struggle then you can move on to one of the other tutorials instead.

## RNA-seq analysis of MRSA
As you might remember from the [intro](index), we are attempting to understand how lytic bacteriophages can be used as a future therapy for the multiresistant bacteria MRSA (methicillin-resistant _Staphylococcus aureus_). In order to do this we have performed RNA-seq of three samples, X treated and Y untreated. We've set up a Snakemake workflow for the RNA-seq analysis and it seems to be running nicely. We'd now like to modify this workflow to make it more flexible and reproducible.

This section will leave a little more up to you compared to the previous one. If you get stuck at some point the final workflow after all the modifications can be found in the `git_jupyter_docker` directory.

Let's start by generating the rule graph so that we get an overview of the workflow.

```bash
snakemake -s Snakefile_mrsa --rulegraph | dot -Tpng > rulegraph_mrsa.png
```

There are two differences in this command compared to the one we've used before. The first is that we're using the `-s` flag to specify which Snakemake workflow to run. We didn't need to do that before since `Snakefile` is the default name. The second is that we don't define a target. In the toy example we used `a_b.txt` as a target, and the wildcards were resolved based on that. How come that we don't need to do that here? It turns out that by default Snakemake targets the first rule in a workflow. By convention we call this rule `all` and let it serve as a rule for aggregating the main outputs of the workflow.

![alt text](rulegraph_mrsa.svg)

Now take some time and look through the workflow file and try to understand how the rules fit together. Use the rule graph as aid. The rules represent a quite standard, although somewhat simplified, workflow for RNA-seq analysis. If you are unfamiliar with the purpose of the different operations (index genome, FastQC and so on), then take a look at the [intro](index).

Also generate the job graph in the same manner. Here you can see that three samples will be downloaded from SRA (Sequence Read Archive); SRR935090, SRR935091, and SRR935092. Those will then be quality controlled with FastQC and aligned to a genome. The QC output will be aggregated with MultiQC and the alignments will be used to generate a count table, i.e a table that show how many reads map to each gene for each sample. This count table is then what the downstream analysis will be based on (in the [Rmarkdown exercise](rmarkdown)).

![alt text](dag_mrsa.svg)

Now try to run the whole workflow. Hopefully you see something like this.

```no-highlight
Building DAG of jobs...
Provided cores: 1
Rules claiming more threads will be scaled down.
Job counts:
	count	jobs
	3	align_to_genome
	1	all
	3	fastqc
	1	generate_count_table
	1	generate_rulegraph
	3	get_SRA_by_accession
	1	get_genome_fasta
	1	get_genome_gff3
	1	index_genome
	1	multiqc
	3	sort_bam
	19

rule get_genome_fasta:
    output: data/raw_external/NCTC8325.fa.gz
    jobid: 18

--2017-11-16 22:13:28--  ftp://ftp.ensemblgenomes.org/pub/bacteria/release-37/fasta/bacteria_18_collection/staphylococcus_aureus_subsp_aureus_nctc_8325/dna//Staphylococcus_aureus_subsp_aureus_nctc_8325.ASM1342v1.dna_rm.toplevel.fa.gz
           => ‘data/raw_external/NCTC8325.fa.gz’
Resolving ftp.ensemblgenomes.org (ftp.ensemblgenomes.org)... 193.62.197.94
Connecting to ftp.ensemblgenomes.org (ftp.ensemblgenomes.org)|193.62.197.94|:21... connected.
Logging in as anonymous ... Logged in!
==> SYST ... done.    ==> PWD ... done.
.
.
[lots of stuff]
.
.
localrule all:
    input: data/final/counts.tsv, results/multiqc.html, results/rulegraph.png
    jobid: 0


Finished job 0.
19 of 19 steps (100%) done
```

After everything is done the workflow will have resulted in a bunch of files in the directories `data`, `intermediate` and `results`. Take some time to look through the structure, in particular the quality control reports in `results` and the count table in `data/final`.

**As noted above, the rest of this tutorial will go over many more nice features of Snakemake. If you got a taste for Snakemake, just continue along. If you want to save this for another day and rather have time to focus on the remaining tutorials, this would be a good point to exit.**

### Logs
As you probably noticed it was difficult to follow how the workflow progressed since some rules printed a lot of output to the terminal. In some cases this also contained important information, such as statistics on the sequence alignments or genome indexing. This could be valuable for example if you later in the project get weird results and want to debug. It's also important from a reproducibility perspective that the "paper trail" describing how the outputs were generated is saved. Luckily, Snakemake has a feature that can help with this. Just as we define `input` and `output` in a rule we can also define `log`.

```python
rule some_rule:
    input:
        "..."
    output:
        "..."
    log:
        "..."
    shell:
        """
        echo 'Converting {input} to {output}' > {log}
        """
```

A log file is not different from any other output file, but it's dealt with a little differently by Snakemake. For example it's shown in the file summary when using `-D`. It's also a good way to clarify the purpose of the file. We probably don't need to save logs for all the rules, only the ones with interesting output.

* `get_genome_fasta` and `get_genome_gff3` would be good to log since they are dependent on downloading files from an external server.
* `multiqc` aggregates quality control data for all the samples into one html report, and the log contains information about which samples were aggregated.
* `index_genome` outputs some statistics about the genome indexing.
* `align_to_genome` outputs important statistics about the alignments. This is probably the most important log to save.

Now add a log file to some or all of the rules above. A good place to save them to would be `results/logs/rule_name/job_name.log`. Be sure to include any wildcards used in the rule in the job name as well, so that you don't end up with identical names for different samples for example.

You also have to specify in the `shell` section of each rule what you want the log to contain. Some of the programs we use send their log information to standard out, some to standard error and some let us specify a log file via a flag. To save some time you can use the info below.

```bash
# Wget has a -o flag for specifying the log file
wget remote_file -O output_file -o {log}

# MultiQC writes to standard error so we redirect with "2>"
multiqc -n output_file input_files 2> {log}

# Bowtie2-build redirects to standard out so we use ">"
bowtie2-build input_file index_dir > {log}

# Bowtie2 writes the main output to standard out and the log info to standard error
bowtie2 -x index_dir -U input_file > output_file 2> {log}
```

Now rerun the whole workflow by using the `-F` flag. Do the logs contain what they should? Note how much easier it it to follow the progression of the workflow when the rules write to logs instead of to the terminal. If you run with `-D` (or `-S` for a simpler version) you will see that the summary table now also contains the log file for each of the files in the workflow.

### Marking files as temporary
It's not uncommon that workflows in bioinformatics contain temporary files that should be kept for some time and then deleted once they are no longer needed. A typical case could be that some operation generates a file, which is then compressed to save space or indexed to make searching faster. There is then no need to save the original output file. Take a look at the job graph for our workflow again. The output from `align_to_genome` is a bam file, which contains information about all the reads for a sample and where they map in the genome. For downstream processing we need this file to be sorted by genome coordinates. This is what the rule `sort_bam` is for. We therefore end up with both `intermediate/{sra_id}.bam` and `intermediate/{sra_id}.sorted.bam`.

In Snakemake we can mark an output file as temporary like this:

```python
output: temp("...")
```

The file will then be deleted as soon as all jobs where it's an input have finished. Now do this for the output of `align_to_genome`. We have to rerun the rule for it to trigger, so use `-R align_to_genome`. It should look something like this:

```no-highlight
.
.
rule sort_bam:
    input: intermediate/SRR935090.bam
    output: intermediate/SRR935090.sorted.bam
    jobid: 2
    wildcards: sra_id=SRR935090

Removing temporary output file intermediate/SRR935090.bam.
Finished job 2.
.
.
```

!!! tip
  Sometimes you may want to trigger removal of temporary files without actually rerunning the jobs. You can then use the `--touch` flag, which will run though the workflow and update the timestamps on all output files without actually executing the code in `shell`.

Snakemake has a number of options for marking files:

* `temp("...")`: The output file should be deleted once it's no longer needed by any rules.
* `protected("...")`: The output file should be write-protected. Typically used to protect files which require a huge amount of computation time from being accidentally deleted.
* `ancient("...")`: The timestamp of the input file is ignored and it's always assumed to be older than any of the output files.
* `touch("...")`: The output file should be "touched", i.e. created or updated, when the rule has finished. Typically used as "flag files" to enforce some rule execution order without real file dependencies.
* `dynamic("{some_wildcard}...")`: This one is a bit tricky. It's used when the number of output files from a rule cannot be determined beforehand. A typical use case could be if you run some clustering analysis and end up with one file per cluster.

### Shadow rules
Take a look at the rule `generate_count_table` below. Since `input.annotation` is compressed it is first unzipped to a temporary file. `htseq-count` then generates a temporary count table, which is finally prepended with a header containing the sample names. Lastly, the two temporary files are deleted.

```python
rule generate_count_table:
    """
    Generate a count table using htseq-count.
    """
    input:
        bams=["intermediate/SRR935090.sorted.bam", "intermediate/SRR935091.sorted.bam", "intermediate/SRR935092.sorted.bam"],
        annotation="data/raw_external/NCTC8325.gff3.gz"
    output:
        "data/final/counts.tsv"
    shell:
        """
        # htseq-count cannot use .gz, so unzip to a temporary file first
        gunzip -c {input.annotation} > tempfile

        # Save the count table as a temporary file and then prepend a header line
        # with the sample names
        htseq-count --format bam --type gene --idattr gene_id {input.bams} tempfile > tempfile2
        echo '{input.bams}' | tr ' ' '\t' | cat - tempfile2 > {output}

        # Remove the temporary files
        rm tempfile tempfile2
        """
```

There are a number of drawbacks with having files that aren't explicitly part of the workflow as input/output files to rules (as the two temporary files here).

* Snakemake cannot clean up these files if the job fails, as it would do for the normal output files.
* If several jobs are run in parallel there is a risk that they write to `tempfile` at the same time. This can lead to very scary results.
* Sometimes we don't know the names of all the files that a program can generate. It is, for example, not unusual that programs leave some kind of error log if something goes wrong.

All of these issues can be dealt with by using the `shadow` option for a rule. Shadow rules result in that each execution of the rule is run in an isolated temporary directory (located in `.snakemake/shadow/`). There are a few option for `shadow`. The most simple is `shadow: "minimal"` (*NOTE: USE SHADOW: "SHALLOW" UNTIL PR ACCEPTED*), which means that the rule is executed in an empty directory where the input files have been symlinked. For the rule below that means that the only file available would be `input.txt`. The shell commands would generate the files `some_other_junk_file` and `output.txt`. Lastly, Snakemake will move the output file (`output.txt`) to its "real" location and remove the whole shadow directory. We therefore never have to think about manually removing `some_other_junk_file`.

```python
rule some_rule:
    input:
        "input.txt"
    output:
        "output.txt"
    shadow: "minimal"
    shell:
        """
        touch some_other_junk_file
        cp {input} {output}
        """
```

Try this out for the rules where we have to "manually" deal with files that aren't tracked by Snakemake (`multiqc`, `index_genome`, `generate_count_table`). Also remove the shell commands that remove temporary files from those rules, as they are no longer needed. Now rerun the workflow and validate that the temporary files don't show up in your working directory.

Some people use the shadow option for almost every rule and some never use it at all. One thing to keep in mind is that it leads to some extra file operations when the outputs are moved to their final location. This is no issue when the `.snakemake` directory is on the same disk as the output directory, but if you're running on a distributed file system and generate many very large files it might be worth considering other options.

*ADD OUTPUT EXAMPLE WHEN ISSUE DEALT WITH*

### Rule targets
So far we have only defined the inputs/outputs of a rule as strings or in some case a list of strings, but Snakemake allows us to be much more flexible than that. Actually, we can use any Python expression or even functions, as long as they return a string or list of strings. Consider the rule `align_to_genome` below.

```python
rule align_to_genome:
    """
    Align a fastq file to a genome index using Bowtie 2.
    """
    input:
        "data/raw_internal/{sra_id}.fastq.gz",
        "intermediate/NCTC8325.1.bt2",
        "intermediate/NCTC8325.2.bt2",
        "intermediate/NCTC8325.3.bt2",
        "intermediate/NCTC8325.4.bt2",
        "intermediate/NCTC8325.rev.1.bt2",
        "intermediate/NCTC8325.rev.2.bt2"
    output:
        "intermediate/{sra_id,\w+}.bam"
    shell:
        """
        bowtie2 -x intermediate/NCTC8325 -U {input[0]} > {output}
        """
```

Here we have seven inputs; the fastq file with the reads and six files with similar file names from the Bowtie 2 genome indexing. We can try to tidy this up by using a Python expression to generate a list of these files instead. If you're familiar with Python you could do this with list comprehensions like this:

```python
input:
    "data/raw_internal/{sra_id}.fastq.gz",
    ["intermediate/NCTC8325.{my_substr}.bt2".format(my_substr=substr) for substr in ["1", "2", "3", "4", "rev.1", "rev.2"]]
```

This will take the elements of the list of substrings one by one, and insert that element in the place of `{my_substring}`. Since this type of aggregating rules are quite common, Snakemake also has a more compact way of achieving the same thing.

```python
input:
    "data/raw_internal/{sra_id}.fastq.gz",
    expand("intermediate/NCTC8325.{my_substr}.bt2", my_substr = ["1", "2", "3", "4", "rev.1", "rev.2"])
```

Now change in the rules `index_genome` and `align_to_genome` to use the `expand()` expression.

In the workflow we decide which samples to run by including the SRR ids in the names of the inputs to the rules `multiqc` and `generate_count_table`. This is a potential source of errors since it's easy to change in one place and forget to change in the other. As we've mentioned before, but not really used so far, Snakemake allows us to use Python code "everywhere". Let's therefore define a list of sample ids and put at the very top of the Snakefile, just before the rule `all`.

```python
SAMPLES = ["SRR935090", "SRR935091", "SRR935092"]
```

Now use `expand()` in `multiqc` and `generate_count_table` to use `SAMPLES` for the sample ids. Much better!

### Parameters
In a typical bioinformatics project considerable efforts are spent on tweaking parameters for the various programs involved. It would be inconvenient if you had to change in the shell scripts themselves everytime you wanted to run with a new setting. Luckily, there is a better option for this; the `params` keyword.

```python
rule some_rule:
    input:
        "..."
    output:
        "..."
    params:
        cutoff=2.5
    shell:
        """
        some_program --cutoff {params.cutoff} {input} {output}
        """
```

In our workflow we run most of the programs with default settings. However, there is one parameter in the rule `get_SRA_by_accession` that we use for determining how many reads we want to retrieve from SRA for each read (`-X 25000`). Change in this rule to use the parameter `max_reads` instead.

```python
rule get_SRA_by_accession:
    """
    Retrieve a single-read FASTQ file from SRA (Sequence Read Archive) by run accession number.
    """
    output:
        "data/raw_internal/{sra_id}.fastq.gz"
    shell:
        """
        fastq-dump {wildcards.sra_id} -X 25000 --readids \
            --dumpbase --skip-technical --gzip -Z > {output}

        # This clears a cache where SRA Tools reserve a lot of space
        cache-mgr --clear >/dev/null 2>&1
        """
```

Now when we have declare `max_reads` as a parameter we can change the value from the command line by using the `snakemake --config [KEY=VALUE [KEY=VALUE ...]]` syntax. Try this out for yourself.

From a reproducibility perspective it's not optimal to set parameters from the command line, since it's difficult to keep track of which parameter values that were used. A much better alternative is to use the `--configfile FILE` option. Here we can collect all the project-specific settings, sample ids, and so on in one file. This also enables us to write the Snakefile in a more general manner so that it can be better reused between projects. Like several other files used in these tutorials, this file should be in [yaml format](https://en.wikipedia.org/wiki/YAML). Create the file below and save it as `config.yml`.

```yaml
max_reads: 25000
```

If we now supply `--configfile config.yml` Snakemake will parse this into a global dictionary called `config` that is available from all rules. We can therefore modify `get_SRA_by_accession` to look something like this. Now try this out yourself.

```python
rule get_SRA_by_accession:
    """
    Retrieve a single-read FASTQ file from SRA (Sequence Read Archive) by run accession number.
    """
    output:
        "data/raw_internal/{sra_id}.fastq.gz"
    params:
        max_reads = config["max_reads"]
    shell:
        """
        fastq-dump {wildcards.sra_id} -X {params.max_reads} --readids \
            --dumpbase --skip-technical --gzip -Z > {output}

        # This clears a cache where SRA Tools reserve a lot of space
        cache-mgr --clear >/dev/null 2>&1
        """
```

!!! tip
  Rather than supplying the config file from the command line you could also add the line `configfile: "config.yml"` to the top of your Snakefile.

If we want to move all project-specific information to `config.yml` and let the Snakefile be a more general RNA-seq analysis workflow we need to:

* Specify which samples to run.
* Specify which genome to align to and where to download its sequence and annotation files.
* (Any other parameters we might need to make it into a general workflow, e.g. to support both paired-end and single-read sequencing)

The first point is straight forward; rather than using `SAMPLES = ["..."]` in the Snakefile we define it as a parameter in `config.yml`. Do a dry-run afterwards to make sure that everything works as expected. You can either add it as a list as it was expressed before, or you can use this yaml notation:

```yaml
sample_ids:
- SRR935090
- SRR935091
- SRR935092
```

The second point is trickier. Writing workflows in Snakemake is quite straightforward when the logic of the workflow is reflected in the file names, i.e. `my_sample.trimmed.deduplicated.sorted.fastq`, but that isn't always the case. In our case we have the FTP paths to the genome sequence and annotation where the naming doesn't quite fit with the rest of the workflow. The easiest solution is probably to make three parameters to hold these values, say `genome_id`, `genome_fasta_path` and `genome_gff_path`, but we will go for a somewhat more complex but very useful alternative. We want to construct a dictionary where something that will be a wildcard in the workflow is the key and the troublesome name is the value. An example might make this clearer (this is also in `config.yml`). This is a nested dictionary where "genomes" is a key with with another dictionary as value, which in turn has genome ids as keys and so on. The idea is that we have a wildcard in the workflow that takes the id of a genome as value (either "NCTC8325" or "ST398" in this case). The fasta and gff3 paths can then be retrieved based on the value of the wildcard.

```yaml
genomes:
  NCTC8325:
    fasta: ftp://ftp.ensemblgenomes.org/pub/bacteria/release-37/fasta/bacteria_18_collection/staphylococcus_aureus_subsp_aureus_nctc_8325/dna//Staphylococcus_aureus_subsp_aureus_nctc_8325.ASM1342v1.dna_rm.toplevel.fa.gz
    gff3: ftp://ftp.ensemblgenomes.org/pub/bacteria/release-37/gff3/bacteria_18_collection/staphylococcus_aureus_subsp_aureus_nctc_8325//Staphylococcus_aureus_subsp_aureus_nctc_8325.ASM1342v1.37.gff3.gz
  ST398:
    fasta: ftp://ftp.ensemblgenomes.org/pub/bacteria/release-37/fasta/bacteria_18_collection//staphylococcus_aureus_subsp_aureus_st398/dna/Staphylococcus_aureus_subsp_aureus_st398.ASM958v1.dna.toplevel.fa.gz
    gff3: ftp://ftp.ensemblgenomes.org/pub/bacteria/release-37/gff3/bacteria_18_collection/staphylococcus_aureus_subsp_aureus_st398//Staphylococcus_aureus_subsp_aureus_st398.ASM958v1.37.gff3.gz
```

Let's now look at how to do the mapping from genome id to fasta path in the rule `get_genome_fasta`. This is how the rule currently looks (if you have added the log section as previously described).

```python
rule get_genome_fasta:
    """
    Retrieve the sequence in fasta format for a genome.
    """
    output:
        "data/raw_external/NCTC8325.fa.gz"
    log:
        "results/logs/get_genome_fasta/NCTC8325.log"
    shell:
        """
        wget ftp://ftp.ensemblgenomes.org/pub/bacteria/release-37/fasta/bacteria_18_collection/staphylococcus_aureus_subsp_aureus_nctc_8325/dna//Staphylococcus_aureus_subsp_aureus_nctc_8325.ASM1342v1.dna_rm.toplevel.fa.gz -O {output} -o {log}
        """
```

We don't want the hardcoded genome id, so replace that with a wildcard, say `{genome_id}`. In general, we want the rules as far downstream as possible in the workflow to be the ones that determine what the wildcards should resolve to. In our case this is `align_to_genome`. You can think of it like the rule that really "needs" the file asks for it, and then it's up to Snakemake to determine how it can use all the available rules to generate it. Here we say "I need this genome index to align my sample to" and it's up to Snakemake to determine how to download and build the index.

We're almost done, but there is one more tricky thing left. Take a look below at the `params` section. Here we have defined a function to generate the parameter `fasta_path`. This is what's called an anonymous function, or lambda expression, but any normal function would work as well. What happens is that the function takes the `wildcards` object as its only argument. This object allows access to the wildcards values via attributes (here `wildcards.genome_id`). The function will then look in the nested `config` dictionary and return the value of the fasta path for the key `wildcards.genome_id`. Neat!

```python
rule get_genome_fasta:
    """
    Retrieve the sequence in fasta format for a genome.
    """
    output:
        "data/raw_external/{genome_id}.fa.gz"
    params:
        fasta_path = lambda wildcards: config["genomes"][wildcards.genome_id]["fasta"]
    log:
        "results/logs/get_genome_fasta/{genome_id}.log"
    shell:
        """
        wget {params.fasta_path} -O {output} -o {log}
        """
```

Now change in `get_genome_gff3` in the same way. Also change in `index_genome` to use a wildcard rather than a hardcoded genome id. Here you will run into a complication if you have followed the previous instructions and use the `expand()` expression. We want the list to expand to `["intermediate/{genome_id}.1.bt2", "intermediate/{genome_id}.2.bt2", ...]`, but for it to do so we must use double curly brackets around the wildcard, i.e. `{{genome_id}}`. Lastly, we need to define somewhere which genome id we actually want to use. This needs to be done both in `align_to_genome` and `generate_count_table`. Do this by introducing a parameter in `config.yml` called "genome_id". See below for an example for `align_to_genome`.

```python
output:
    expand("intermediate/{genome_id}.{substr}.bt2", genome_id = config["genome_id"], substr = ["1", "2", "3", "4", "rev.1", "rev.2"])
```

### Summary
Well done!

* You have a general RNA-seq pipeline which can easily be reused between projects, thanks to clear separation between code and settings.
* You have great traceability due to logs and summary tables.
* You have clearly defined the environment for the workflow using Conda.
* You have kept the workflow neat and free from temporary files by using `temp()` and `shadow`.
* You have a logical directory structure which makes it easy to separate raw data, intermediate files, and results.
* You have set up you project in a way that makes it very easy to distribute and reproduce, either via git, via Snakemake's `--archive` option, or via a Docker image.

# Where to go next?

# TODO
run/script
list-x-changes (mention version)
cluster
conda/docker
Archive
