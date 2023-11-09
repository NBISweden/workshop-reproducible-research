Here are some useful resources if you want to read more about Jupyter in
general:

* The [Jupyter project site](http://jupyter.org) contains a lot of information
  and inspiration.
* The [Jupyter Notebook documentation](
  https://jupyter-notebook.readthedocs.io/en/stable/).
* A [guide](http://ipywidgets.readthedocs.io/en/stable/index.html) to using
  widgets for creating interactive notebooks.

## Running Jupyter notebooks on a cluster

* Login to Uppmax, making sure to use a specific login node, _e.g._ `rackham1`:

```
ssh <your-user-name>@rackham1.uppmax.uu.se
```

* Create/activate a Conda environment containing `jupyter`, _e.g._:

```
mamba create -n jupyter -c conda-forge jupyter
```

* activate the environment, then run:

```
jupyter notebook
```

When the Jupyter server starts up you should see something resembling:
```
[I 11:00:00.000 NotebookApp] Serving notebooks from local directory: <path-to-your-local-dir>
[I 11:00:00.000 NotebookApp] Jupyter Notebook 6.4.6 is running at:
[I 11:00:00.000 NotebookApp] http://localhost:8889/?token=357d65100058efa40a0641fce7005addcff339876c5e8000
[I 11:00:00.000 NotebookApp]  or http://127.0.0.1:8889/?token=357d65100058efa40a0641fce7005addcff339876c5e8000
[I 11:00:00.000 NotebookApp] Use Control-C to stop this server and shut down all kernels (twice to skip confirmation).
```

Now a Jupyter notebook server is running on the Uppmax end. The line that says:
```
[I 11:00:00.000 NotebookApp] http://localhost:8889/?token=357d65100058efa40a0641fce7005addcff339876c5e8000
```

Contains information on the port used on the server side (8889 in this case) and
the token required to use the server (`357d65100058efa40a0641fce7005addcff339876c5e8000`).

Next step is to use this information to login to the server from your local
computer.

**On your local computer**

In a terminal, run the following command to start port forwarding of
port 8080 on your local computer to the remote port on the Uppmax side. Replace
<remote-port> with the port given when you started the server on Uppmax. Also
replace <your-user-name> with your user name on Uppmax.

```bash
ssh -N -L localhost:8080:localhost:<remote-port> <your-user-name>@rackham1.uppmax.uu.se
```

As long as this process is running the port forwarding is running. To disable it
simply interrupt it with `CTRL + C`.

Connect to the Jupyter server by opening `localhost:8080` in your browser. When
prompted, paste the token you got when starting the server on Uppmax.

You are now (hopefully) accessing the Jupyter server that's running on Uppmax,
via your local browser.

## Integrating notebooks with Snakemake workflows

In the [case study](jupyter-8-the-mrsa-case-study) section of this tutorial we
created a Jupyter notebook that used output from a Snakemake workflow
and produced some summary results and plots. Wouldn't it be nice if this was
actually part of the workflow itself? To generate a HTML version of the notebook
we can use what we learned in the section about
[converting notebooks](jupyter-7-converting-notebooks). The command to execute the notebook
and save it in HTML format in a file `results/supplementary.html` would be:

```bash
jupyter nbconvert --to HTML --output-dir results --output supplementary.html --execute supplementary_material.ipynb
```

This command could be used in a rule, _e.g._ `make_supplementary`, the input of
which would be `results/tables/counts.tsv`, `intermediate/multiqc_general_stats.txt`,
and `results/rulegraph.png`. See if you can work out how to implement such a
rule at the end of the `Snakefile` found in the `jupyter/` directory. You can
find an example in the code chunk below:

```python
rule make_supplementary:
    output:
            "results/supplementary.html"
    input:
        counts = "results/tables/counts.tsv",
        summary = "results/tables/counts.tsv.summary",
        multiqc_file = "intermediate/multiqc_general_stats.txt",
        rulegraph = "results/rulegraph.png"
    params:
        base = lambda wildcards, output: os.path.basename(output[0]),
        dir = lambda wildcards, output: os.path.dirname(output[0])
    shell:
        """
        jupyter nbconvert --to HTML --output-dir {params.dir} --output {params.base} \
            --execute supplementary_material.ipynb
        """
```

> **Note** <br>
> The Conda environment for the Jupyter tutorial does not contain the
> Snakemake package so if you wish to test the rule _e.g._ by running
> `snakemake -c 1 results/supplementary.html` you first install Snakemake
> into the environment with `mamba install snakemake`.

## More integrations

Snakemake actually supports the execution of notebooks via the `notebook:`
rules directive. See more about Jupyter integration in the
[Snakemake docs](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#jupyter-notebook-integration).
In the `notebook:` directive of such a rule you specify the path to a Jupyter
notebook (relative to the Snakefile) which is then executed when
the rule is run.

So how is this useful?

In the notebook itself this gives you access to a `snakemake` object
containing information about **input** and **output** files for the rule via
`snakemake.input` and `snakemake.output`. Similarly, you can access rule
**wildcards** with `snakemake.wildcards`, **params** with `snakemake.params`,
and **config** settings with `snakemake.config`.

When Snakemake runs the rule with the `notebook:` directive `jupyter-nbconvert`
is used to execute the notebook. No HTML output is generated here but it is
possible to store a version of the notebook in its final processed form by
adding the following to the rule:

```python
log:
    notebook="<path>/<to>/<processed>/<notebook.ipynb>"
```

Because you won't get the notebook in full HTML glory, this type of
integration is better suited if you want to use a notebook to generate figures
and store these in local files (_e.g._ pdf/svg/png formats).

We'll use the `supplementary_material.ipynb` notebook as an example! Let's say
that instead of exporting the entire notebook to HTML we want a rule that
outputs pdf versions of the barplot and heatmap figures we created.

Let's start by setting up the rule. For simplicity we'll use the same input as
when we edited the notebook in the first place. The output will be
`results/barplot.pdf` and `results/heatmap.pdf`. Let's also output a finalized
version of the notebook using the `log: notebook=` directive:

```python
rule make_supplementary_plots:
    output:
        barplot = "results/barplot.pdf",
        heatmap = "results/heatmap.pdf"
    input:
        counts="results/tables/counts.tsv",
        summary="results/tables/counts.tsv.summary",
        multiqc="intermediate/multiqc_general_stats.txt",
        rulegraph="results/rulegraph.png"
    log:
        notebook = "results/supplementary.ipynb"
```

The notebook will now have access to `snakemake.input.counts`,
`snakemake.output.barplot` and `snakemake.output.heatmap` when executed from
within the workflow. Let's go ahead and edit the notebook! In the cell where we
defined notebook parameters edit the code so that it looks like this:

```python
counts_file=snakemake.input.counts
summary_file=snakemake.input.summary
multiqc_file=snakemake.input.multiqc
rulegraph_file=snakemake.input.rulegraph

SRR_IDs=snakemake.params.SRR_IDs
GSM_IDs=snakemake.params.GSM_IDs
GEO_ID=snakemake.params.GEO_ID
```

Notice that we set the `SRR_IDs`, `GSM_IDs` and `GEO_ID` variables using
variables in `snakemake.params`? However, we haven't defined these in our rule
yet so let's go ahead and do that now. Add the `params` section so that the
`make_supplementary_plots` in the Snakefile looks like this:

```python
rule make_supplementary_plots:
    output:
        barplot = "results/barplot.pdf",
        heatmap = "results/heatmap.pdf"
    input:
        counts="results/tables/counts.tsv",
        summary="results/tables/counts.tsv.summary",
        multiqc="intermediate/multiqc_general_stats.txt",
        rulegraph="results/rulegraph.png"
    log:
        notebook = "results/supplementary.ipynb"
    params:
        SRR_IDs = ["SRR935090","SRR935091","SRR935092"],
        GSM_IDs = ["GSM1186459", "GSM1186460", "GSM1186461"],
        GEO_ID = "GSE48896"
    notebook: "supplementary_material.ipynb"
```

> **Tip** <br>
> One way to further generalize this rule could be to define the SRR_IDs, GSM_IDs
> and GEO_ID parameters in a config file instead, in which case they would be
> directly accessible from within the notebook using `snakemake.config['SRR_IDs']`
> etc.

Now the rule contains everything needed, but we still need to edit the notebook
to save the plots to the output files. First, edit the cell that generates the
barplot so that it looks like this:

```python
# Create a stacked barplot
ax = summary_plot_data.T.plot(kind="bar", stacked=True, color=colors)
# Move legend and set legend title
ax.legend(bbox_to_anchor=(1,1), title="Category");
plt.savefig(snakemake.output.barplot, dpi=300, bbox_inches="tight") ## <-- Add this line!
```

Finally, edit the cell that generates the heatmap so that it looks like this:

```python
count_data = counts.loc[:, SRR_IDs]
heatmap_data = count_data.loc[(count_data.std(axis=1).div(count_data.mean(axis=1))>1.2)&(count_data.max(axis=1)>5)]
heatmap_data = heatmap_data.rename(columns = name_dict['title'])
with sns.plotting_context("notebook", font_scale=0.7):
    ax = sns.clustermap(data=np.log10(heatmap_data+1), cmap="RdYlBu_r",
                        method="complete", yticklabels=True, linewidth=.5,
                        cbar_pos=(.7, .85, .05, .1), figsize=(3,9))
    plt.setp(ax.ax_heatmap.get_xticklabels(), rotation=270)
    plt.savefig(snakemake.output.heatmap, dpi=300, bbox_inches="tight") ## <-- Add this line!
```

Now you can run the following to generate the plots:

```bash
snakemake -c 1 make_supplementary_plots
```

## Presentations with Jupyter

As if all the above wasn't enough you can also create presentations/slideshows
with Jupyter! Simply use mamba to install the
[RISE](https://rise.readthedocs.io/en/stable/) extension to your Jupyter
environment:

```bash
mamba install -c conda-forge rise
```

> **Attention!** <br>
> Unfortunately, the RISE extension is currently not supported in Jupyter
> lab so you are forced to use the classic notebook interface here.

Then open up a notebook of your choice. In the menu click **View** -> **Cell
Toolbar** -> **Slideshow**. Now every cell will have a drop-down in the upper
right corner allowing you to set the cell type:

- **Slide**: a regular slide
- **Sub-Slide**: a regular slide that will be displayed _below_ the previous
- **Fragment**: these cells split up slides so that content (fragments) are
  added only when you press Space
- **Skip**: these cells will not appear in the presentation
- **Notes**: these cells act as notes, shown in the speaker view but not in the
  main view

The presentation can be run directly from the notebook by clicking the
'Enter/Exit RISE Slideshow' button (looks like a bar chart) in the toolbar, or
by using the keyboard shortcut `Alt-r`. Running it directly from a notebook
means you can also edit and run cells during your presentation. The downside is
that the presentation is not as portable because it may rely on certain software
packages that other people are not comfortable with installing.

You can also export the notebook to an HTML-file with `jupyter nbconvert
--execute --to SLIDES <your-notebook.ipynb>`. The resulting file, with the
slideshow functionality included, can be opened in any browser. However, in
this format you cannot run/edit cells.

## Parameterising notebooks

- Use [papermill](https://papermill.readthedocs.io/en/latest/) to parameterise
  notebooks and run them as scripts.