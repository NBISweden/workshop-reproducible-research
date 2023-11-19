You may have noticed that a lot of the functionality in Jupyter is overlapping
with Quarto. And you may be wondering which one to use. This is a difficult
question to answer as it will depend on your use-case and personal preference.
As such, any answer will be subjective, but we'll try to give you some 
pointers on how to get the best out of both worlds.

While similar in some ways Jupyter and Quarto are not completely overlapping.
Quarto is great for generating high-quality reports and manuscripts, and is
agnostic to the programming language used. Jupyter on the other hand is
great for interactive data analysis and exploration with a more direct
connection between code and output. While Jupyter is also somewhat agnostic to
programming language, it is most commonly used with Python and with both the
Jupyter and Python ecosystem at its back it can be customized with a lot of
different extensions and plugins.

The good news is that the two can be used together allowing you to get the best
of both. For example, you may like the professional look of rendered Quarto
documents but really like the interactive and exploratory nature of Jupyter.
Well you can simply work as you normally do in Jupyter and then use Quarto to
render the notebook to a high-quality report or manuscript.

To give you an example, take a look at the `supplementary_material.ipynb` file
in the `jupyter/` tutorial directory. Open this notebook in the Jupyter lab
interface (make sure you have activated the `jupyter-env` Conda environment).

As you can see this notebook contains some brief descriptions in Markdown and
code to generate a few plots. It uses the output from the MRSA case-study
Snakemake workflow you worked on in the Snakemake tutorial. This is a common
use-case for Jupyter notebooks; to generate summary statistics and plots from
the results of a workflow run. (A real-world example could of course include a
lot more in-depth exploratory analyses).

Now, let's say you want to share the results of this notebook with your PI or
collaborators. We could simply share the notebook file, or as we saw in the
previous section, convert it to HTML or PDF via `jupybter nbconvert`. 

Let's do that first so we have something to compare with. Run the following:

```bash
jupyter nbconvert --to HTML --output supplementary_material.nbconvert.html supplementary_material.ipynb
```

Open the `supplementary_material.nbconvert.html` file in a browser to see that
it looks like you expect. This looks more or less like the original notebook.

Now let's go one step further and render the notebook to a high-quality report
using Quarto. We can actually add a YAML header to the notebook with some
document options that Quarto understands. Create a new cell in the notebook
(from the Jupyter lab interface) and move it to the top. In this cell, add the
following:

```yaml
---
title: Supplementary material
subtitle: Supplementary tables and plots for the MRSA study
format:
    html:
        embed-resources: true
        code-fold: true
        code-tools: true
language:
  code-summary: Click to show code 
bibliography: references.bib
---
```

Set the cell type to `Markdown`, then run the cell. Most likely that cell will
look rather weird but that's OK. We'll fix that in a bit.

Save the notebook and now render the document with Quarto from the commandline:

```bash
quarto render supplementary_material.ipynb
```

Open up the `supplementary_material.html` file in a browser and compare it to 
the `supplementary_material.nbconvert.html` file. You should see that the
Quarto version looks a lot better. The fact that Quarto supports rendering of
Jupyter notebooks means you can keep editing your notebooks as you normally
would and use Quarto for rendering the final document. Also there's very little
we had to change in the notebook to make it work with Quarto. If you look
closely at the code cells used to generate the plots and table you'll see that
they contain code-chunk options in the same form we used in the Quarto tutorial.
These options do not impact the notebook when run in Jupyter, making it easy to
use the two tools in combination.

Let's go back to the YAML header cell and fix how it looks in the Jupyter
notebook. The reason it looks weird is that Jupyter doesn't understand the
syntax. But luckily there's a Jupyter lab Quarto extension you can install to
fix this. Click the extension icon in the left sidebar and search for `quarto`.
Install the `jupyterlab-quarto` extension and then reload the page. Now the YAML
header should look a lot better.

Try adding more options to the header to customize the look of the rendered
document. For instance you could:

- add a Table of contents with (`toc: true`)
- try out different
  [themes](https://quarto.org/docs/output-formats/html-themes.html)
- add your name as author (`author: Your Name`)
- add a date (`date: last-modified`)

and much more.

> **Quick recap** <br>
> In this section we've learned:
>
> - How to use Quarto to render Jupyter notebooks to high-quality reports.
