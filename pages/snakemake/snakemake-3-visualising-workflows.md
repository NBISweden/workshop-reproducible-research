All that we've done so far could quite easily be done in a simple shell script
that takes the input files as parameters. Let's now take a look at some of the
features where a WfMS like Snakemake really adds value compared to a more
straightforward approach. One such feature is the possibility to visualize your
workflow. Snakemake can generate three types of graphs, one that shows how the
rules are connected, one that shows how the jobs (*i.e.* an execution of
a rule with some given inputs/outputs/settings) are connected, and finally one
that shows rules with their respective input/output files.

First we look at the rule graph. The following command will generate a rule graph
in the dot language and pipe it to the program `dot`, which in turn will save
a visualization of the graph as a PNG file (if you're having troubles displaying
PNG files you could use SVG or JPG instead).

```bash
snakemake --rulegraph a_b.txt | dot -Tpng > rulegraph.png
```

![](images/rulegraph.svg)

This looks simple enough, the output from the rule `convert_to_upper_case` will
be used as input to the rule `concatenate_files`.

For a more typical bioinformatics project it can look something like this when you
include all the rules from processing of the raw data to generating figures for the paper.

![](images/rulegraph_complex.svg)

While saying that it's easy to read might be a bit of a stretch, it definitely
gives you a better overview of the project than you would have without a WfMS.

The second type of graph is based on the jobs, and looks like this for our
little workflow (use `--dag` instead of `--rulegraph`).

```bash
snakemake --dag a_b.txt | dot -Tpng > jobgraph.png
```

![](images/jobgraph.svg)

The main difference here is that now each node is a job instead of a rule. You
can see that the wildcards used in each job are also displayed. Another
difference is the dotted lines around the nodes. A dotted line is Snakemake's
way of indicating that this rule doesn't need to be rerun in order to generate
`a_b.txt`. Validate this by running `snakemake -n -r a_b.txt` and it should say
that there is nothing to be done.


We've discussed before that one of the main purposes of using a WfMS is that it
automatically makes sure that everything is up to date. This is done by
recursively checking that outputs are always newer than inputs for all the
rules involved in the generation of your target files. Now try to change the
contents of `a.txt` to some other text and save it. What do you think will
happen if you run `snakemake -n -r a_b.txt` again?

??? example "Click to show the solution"
    ```no-highlight
    $ snakemake -n -r a_b.txt

    Building DAG of jobs...
    Job stats:
    job                      count    min threads    max threads
    ---------------------  -------  -------------  -------------
    concatenate_files            1              1              1
    convert_to_upper_case        1              1              1
    total                        2              1              1


    [Mon Oct 25 17:00:02 2021]
    rule convert_to_upper_case:
        input: a.txt
        output: a.upper.txt
        jobid: 1
        reason: Updated input files: a.txt
        wildcards: some_name=a
        resources: tmpdir=/var/folders/p0/6z00kpv16qbf_bt52y4zz2kc0000gp/T


    [Mon Oct 25 17:00:02 2021]
    rule concatenate_files:
        input: a.upper.txt, b.upper.txt
        output: a_b.txt
        jobid: 0
        reason: Input files updated by another job: a.upper.txt
        wildcards: first=a, second=b
        resources: tmpdir=/var/folders/p0/6z00kpv16qbf_bt52y4zz2kc0000gp/T

    Job stats:
    job                      count    min threads    max threads
    ---------------------  -------  -------------  -------------
    concatenate_files            1              1              1
    convert_to_upper_case        1              1              1
    total                        2              1              1

    This was a dry-run (flag -n). The order of jobs does not reflect the order of execution.
    ```

Were you correct? Also generate the job graph and compare to the one generated
above. What's the difference? Now rerun without `-n` and validate that
`a_b.txt` contains the new text (don't forget to specify `-c 1`). Note that
Snakemake doesn't look at the contents of files when trying to determine what has
changed, only at the timestamp for when they were last modified.

We've seen that Snakemake keeps track of if files in the workflow have changed,
and automatically makes sure that any results depending on such files are
regenerated. What about if the rules themselves are changed? It turns out that
since version 7.8.0 Snakemake keeps track of this automatically.

Let's say that we want to modify the rule `concatenate_files` to also include
which files were concatenated.

```python
rule concatenate_files:
    output:
        "{first}_{second}.txt"
    input:
        "{first}.upper.txt",
        "{second}.upper.txt"
    shell:
        """
        echo 'Concatenating {input}' | cat - {input[0]} {input[1]} > {output}
        """
```

!!! Note
    It's not really important for the tutorial, but the shell command used here
    first outputs "Concatenating " followed by a space delimited list of the
    files in `input`. This string is then sent to the program `cat` where it's
    concatenated with `input[0]` and `input[1]` (the parameter `-` means that
    it should read from standard input). Lastly, the output from `cat` is sent
    to `{output}`.

If you now run the workflow as before you should see:
```bash
rule concatenate_files:
    input: a.upper.txt, b.upper.txt
    output: a_b.txt
    jobid: 0
    reason: Code has changed since last execution
    wildcards: first=a, second=b
```

because although no files involved in the workflow have been changed, Snakemake
recognizes that the workflow code itself has been modified and this triggers
a re-run.

Snakemake is aware of changes to four categories of such "rerun-triggers":
"input" (changes to rule input files), "params" (changes to the rule `params` section),
"software-env" (changes to conda environment files specified by the `conda:`
directive) and "code" (changes to code in the `shell:`, `run:`, `script:` and
`notebook:` directives).

Prior to version 7.8.0, only changes to the modification time of input files would
trigger automatic re-runs. To run Snakemake with this previous behaviour you
can use the setting `--rerun-triggers mtime` at the command line.
Change the `shell:` section of the `concatenate_files` rule back to the previous
version, then try running: `snakemake -n -r a_b.txt --rerun-triggers mtime` and
you should again see `Nothing to be done (all requested files are present and up to date).`

You can also export information on how all files were generated (when, by which
rule, which version of the rule, and by which commands) to a tab-delimited file
like this:

```bash
snakemake a_b.txt -c 1 -D > summary.tsv
```

The content of `summary.tsv` is shown in the table below:

<table class="table table-hover table-condensed" border=1; style="margin-left:auto; margin-right:auto;">
    <thead style="background-color:{{config.extra.color_table_header}}">
        <tr>
            <td style="padding:5px"> <font size="3"><b> output_file </b> </td>
            <td style="padding:5px"> <font size="3"><b> date </b> </td>
            <td style="padding:5px"> <font size="3"><b> rule </b> </td>
            <td style="padding:5px"> <font size="3"><b> version </b> </td>
            <td style="padding:5px"> <font size="3"><b> log-file(s) </b> </td>
            <td style="padding:5px"> <font size="3"><b> input-file(s) </b> </td>
            <td style="padding:5px"> <font size="3"><b> shellcmd </b> </td>
            <td style="padding:5px"> <font size="3"><b> status </b> </td>
            <td style="padding:5px"> <font size="3"><b> plan </b> </td>
        </tr>
    </thead>
    <tr>
        <td style="padding:5px"> <font size="3"> a_b.txt </td>
        <td style="padding:5px"> <font size="3"> Mon Oct 25 17:01:46 2021 </td>
        <td style="padding:5px"> <font size="3"> concatenate_files </td>
        <td style="padding:5px"> <font size="3"> - </td>
        <td style="padding:5px"> <font size="3"> </td>
        <td style="padding:5px"> <font size="3"> a.upper.txt,b.upper.txt </td>
        <td style="padding:5px"> <font size="3"> cat a.upper.txt b.upper.txt > a_b.txt </td>
        <td style="padding:5px"> <font size="3"> rule implementation changed </td>
        <td style="padding:5px"> <font size="3"> update pending </td>
    </tr>
    <tr>
        <td style="padding:5px"> <font size="3"> a.upper.txt</td>
        <td style="padding:5px"> <font size="3"> Mon Oct 25 17:01:46 2021 </td>
        <td style="padding:5px"> <font size="3"> convert_to_upper_case </td>
        <td style="padding:5px"> <font size="3"> - </td>
        <td style="padding:5px"> <font size="3"> </td>
        <td style="padding:5px"> <font size="3"> a.txt </td>
        <td style="padding:5px"> <font size="3"> tr [a-z] [A-Z] < a.txt > a.upper.txt </td>
        <td style="padding:5px"> <font size="3"> ok </td>
        <td style="padding:5px"> <font size="3"> no update  </td>
    </tr>
    <tr>
        <td style="padding:5px"> <font size="3"> b.upper.txt</td>
        <td style="padding:5px"> <font size="3"> Mon Oct 25 17:01:46 2021 </td>
        <td style="padding:5px"> <font size="3"> convert_to_upper_case </td>
        <td style="padding:5px"> <font size="3"> - </td>
        <td style="padding:5px"> <font size="3"> </td>
        <td style="padding:5px"> <font size="3"> b.txt </td>
        <td style="padding:5px"> <font size="3"> tr [a-z] [A-Z] < b.txt > b.upper.txt </td>
        <td style="padding:5px"> <font size="3"> ok </td>
        <td style="padding:5px"> <font size="3"> no update  </td>
    </tr>
</table>

You can see in the second last column that the rule implementation for `a_b.txt`
has changed. The last column shows if Snakemake plans to regenerate the files
when it's next executed. You can see that for the `concatenate_files` the plan
is `update pending` because we generated the summary with the default behaviour
of using all rerun-triggers.

You might wonder where Snakemake keeps track of all these things? It stores all
information in a hidden subdirectory called `.snakemake`. This is convenient
since it's easy to delete if you don't need it anymore and everything is
contained in the project directory. Just be sure to add it to `.gitignore` so
that you don't end up tracking it with git.

By now you should be familiar with the basic functionality of Snakemake, and
you can build advanced workflows with only the features we have discussed here.
There's a lot we haven't covered though, in particular when it comes to making
your workflow more reusable. In the following section we will start with
a workflow that is fully functional but not very flexible. We will then
gradually improve it, and at the same time showcase some Snakemake features
we haven't discussed yet. Note that this can get a little complex at times, so
if you felt that this section was a struggle then you could move on to one of
the other tutorials instead.

!!! Success "Quick recap"
    In this section we've learned:

    - How to use `--dag` and `--rulegraph` for visualizing the job and rule
      graphs, respectively.
    - How Snakemake reruns relevant parts of the workflow after
      there have been changes.
    - How Snakemake tracks changes to files and code in a workflow
