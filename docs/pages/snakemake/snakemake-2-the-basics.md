In this part of the tutorial we will create a very simple workflow from
scratch, in order to show the fundamentals of how Snakemake works. The workflow
will take two files as inputs, `a.txt` and `b.txt`, and the purpose is to
convert the text in the files to upper case and then to concatenate them.

Run the following shell commands. The first one will make an empty file named
`Snakefile`, which will later contain the workflow. The second and third
commands generate two files containing some arbitrary text.

```bash
touch Snakefile
echo "This is a.txt" > a.txt
echo "This is b.txt" > b.txt
```

Then open `Snakefile` in your favorite text editor. A Snakemake workflow is based on
rules which take some file(s) as input, performs some type of operation on
them, and generate some file(s) as outputs. Here is a very simple rule that
produces `a.upper.txt` as an output, using `a.txt` as input. Copy this rule to
your `Snakefile` and save it.

```python
rule convert_to_upper_case:
    output:
        "a.upper.txt"
    input:
        "a.txt"
    shell:
        """
        tr [a-z] [A-Z] < {input} > {output}
        """
```

!!! warning
    Indentation is important in Snakefiles, so make sure that you have the
    correct number of spaces before `input`/`output`/`shell` and their
    respective subsections. The number of spaces per level doesn't matter as
    long as you're consistent. Here we use four, but you could just as well use
    two for a more compact look. Don't use tabs (unless your editor
    automatically converts them to spaces).

Rules can be given names, here it's `convert_to_upper_case`. While rule names
are not strictly necessary we encourage you to use them and to make an effort to
name your rules in a way that makes it easy to understand the purpose of the rule,
as rule names are one of the main ways to interact with the workflow. The
`shell` section (or directive) contains the shell commands that will convert the
text in the input file to upper case and send it to the output file. In the shell
command string, we can refer to elements of the rule via curly brackets. Here, we
refer to the output file by specifying `{output}` and to the input file by
specifying `{input}`. If you're not very familiar with Bash, this particular
command can be read like "send the contents of `a.txt` to the program `tr`, which
will convert all characters in the set `[a-z]` to the corresponding character in
the set `[A-Z]`, and then send the output to `a.upper.txt`".

Now let's run our first Snakemake workflow. When a workflow is executed
Snakemake tries to generate a set of target files. Target files can be
specified via the command line (or, as you will see later, in several other
ways). Here we ask Snakemake to make the file `a.upper.txt`. It's good practice
to first run with the flag `-n` (or `--dry-run`), which will show what Snakemake
plans to do without actually running anything, and you also need to specify
how many cores to be used for the workflow with `--cores` or `-c`. For now, you
only need 1 so set `-c 1`. You can also use the flag `-p`, for showing the
shell commands that it will execute, and the flag `-r` for showing the reason
for running a specific rule. `snakemake --help` will show you all available
flags.

```no-highlight
$ snakemake -n -c 1 -r -p a.upper.txt

Building DAG of jobs...
Job stats:
job                      count    min threads    max threads
---------------------  -------  -------------  -------------
convert_to_upper_case        1              1              1
total                        1              1              1


[Mon Oct 25 16:48:43 2021]
rule convert_to_upper_case:
    input: a.txt
    output: a.upper.txt
    jobid: 0
    reason: Missing output files: a.upper.txt
    resources: tmpdir=/var/folders/p0/6z00kpv16qbf_bt52y4zz2kc0000gp/T


        tr [a-z] [A-Z] < a.txt > a.upper.txt

Job stats:
job                      count    min threads    max threads
---------------------  -------  -------------  -------------
convert_to_upper_case        1              1              1
total                        1              1              1

This was a dry-run (flag -n). The order of jobs does not reflect the order of execution.
```

You can see that Snakemake plans to run one job: the rule `convert_to_upper_case`
with `a.txt` as input and `a.upper.txt` as output. The reason for doing this is
that it's missing the file `a.upper.txt`. Now execute the workflow without the
`-n` flag and check that the contents of `a.upper.txt` is as expected. Then try
running the same command again. What do you see? It turns out that Snakemake
only reruns jobs if there have been changes to either **the input files, or the
workflow itself**. This is how Snakemake ensures that everything in the
workflow is up to date. We will get back to this shortly.

What if we ask Snakemake to generate the file `b.upper.txt`?

```no-highlight
$ snakemake -n -c 1 -r -p b.upper.txt

Building DAG of jobs...
MissingRuleException:
No rule to produce b.upper.txt (if you use input functions make sure that they don't raise unexpected exceptions).
```

That didn't work well. We could copy the rule to make a similar one for
`b.txt`, but that would be a bit cumbersome. Here is where named wildcards come
in; one of the most powerful features of Snakemake. Simply change the input
from `input: "a.txt"` to `input: "{some_name}.txt"` and the output to `output:
"{some_name}.upper.txt"`. Now try asking for `b.upper.txt` again.

Tada! What happens here is that Snakemake looks at all the rules it has
available (actually only one in this case) and tries to assign values to all
wildcards so that the targeted files can be generated. In this case it was
quite simple, you can see that it says that `wildcards: some_name=b`, but for
large workflows and multiple wildcards it can get much more complex. Named
wildcards is what enables a workflow (or single rules) to be efficiently
generalized and reused between projects or shared between people.

It seems we have the first part of our workflow working, now it's time to make
the second rule for concatenating the outputs from `convert_to_upper_case`. The
rule structure will be similar; the only difference is that here we have two
inputs instead of one. This can be expressed in two ways, either with named
inputs like this:

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

!!! Warning
    If you have multiple inputs or outputs they need to be delimited with
    a comma (as seen above). This is a very common mistake when writing
    Snakemake workflows. The parser will complain, but sometimes the error
    message can be difficult to interpret.

Now try to construct this rule yourself and name it `concatenate_a_and_b`.
The syntax for concatenating two files in Bash is
`cat first_file.txt second_file.txt > output_file.txt`. Call the output `c.txt`.
Run the workflow in Snakemake and validate that the output looks as expected.

Wouldn't it be nice if our workflow could be used for _any_ files, not just
`a.txt` and `b.txt`? We can achieve this by using named wildcards (or in other
ways as we will discuss later). As we've mentioned, Snakemake looks at all the
rules it has available and tries to assign values to all wildcards so that the
targeted files can be generated. We therefore have to name the output file in
a way so that it also contains information about which input files it should be
based on. Try to figure out how to do this yourself. If you're stuck you can
look at the spoiler below, but spend some time on it before you look. Also
rename the rule to `concatenate_files` to reflect its new more general use.

??? example "Click to show the solution"
    ```python
    rule concatenate_files:
        output:
            "{first}_{second}.txt"    
        input:
            "{first}.upper.txt",
            "{second}.upper.txt"
        shell:
            """
            cat {input[0]} {input[1]} > {output}
            """
    ```

We can now control which input files to use by the name of the file we ask
Snakemake to generate. Run the workflow without the flag `-n` (or `--dry-run`)
to execute both rules, providing one core with `-c 1` (or `--cores 1`):

```no-highlight
$ snakemake a_b.txt -c 1

Building DAG of jobs...
Using shell: /bin/bash
Provided cores: 1 (use --cores to define parallelism)
Rules claiming more threads will be scaled down.
Job stats:
job                      count    min threads    max threads
---------------------  -------  -------------  -------------
concatenate_files            1              1              1
convert_to_upper_case        2              1              1
total                        3              1              1

Select jobs to execute...

[Mon Oct 25 16:51:52 2021]
rule convert_to_upper_case:
    input: b.txt
    output: b.upper.txt
    jobid: 2
    wildcards: some_name=b
    resources: tmpdir=/var/folders/p0/6z00kpv16qbf_bt52y4zz2kc0000gp/T

[Mon Oct 25 16:51:53 2021]
Finished job 2.
1 of 3 steps (33%) done
Select jobs to execute...

[Mon Oct 25 16:51:53 2021]
rule convert_to_upper_case:
    input: a.txt
    output: a.upper.txt
    jobid: 1
    wildcards: some_name=a
    resources: tmpdir=/var/folders/p0/6z00kpv16qbf_bt52y4zz2kc0000gp/T

[Mon Oct 25 16:51:53 2021]
Finished job 1.
2 of 3 steps (67%) done
Select jobs to execute...

[Mon Oct 25 16:51:53 2021]
rule concatenate_files:
    input: a.upper.txt, b.upper.txt
    output: a_b.txt
    jobid: 0
    wildcards: first=a, second=b
    resources: tmpdir=/var/folders/p0/6z00kpv16qbf_bt52y4zz2kc0000gp/T

[Mon Oct 25 16:51:53 2021]
Finished job 0.
3 of 3 steps (100%) done
```

Neat!

!!! Tip
    You can name a file whatever you want in a Snakemake workflow, but you will
    find that everything falls into place much nicer if the filename reflects
    the file's path through the workflow, *e.g.* `sample_a.trimmed.deduplicated.sorted.bam`.

!!! Success "Quick recap"
    In this section we've learned:

    - How a simple Snakemake rule looks.
    - How to define target files when executing a workflow.
    - How to use named wildcards for writing generic and flexible rules.
