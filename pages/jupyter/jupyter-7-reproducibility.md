Now that you have a feeling for what Jupyter can do we'll spend a
little time on things to consider specifically from a repdroducibility point of
view when it comes to Jupyter notebooks.


## Version control of Jupyter notebooks

As we've seen, Jupyter notebooks are plain-text JSON files. This means that they
can be version controlled with Git just like any other text file. However,
because of the way Jupyter notebooks store their content, the diffs produced by
Git can be difficult to interpret. Luckily, there are tools that can provide
content-aware diffs and merge functionality for Jupyter notebooks.

One such tool is [nbdime](https://nbdime.readthedocs.io/en/latest/). `nbdime` is
built to understand the structure of Jupyter notebooks and can therefore
generate diffs that are easier to read. It can also be used to merge changes
made to notebooks, which is great especially when collaborating on notebooks
with others.

`nbdime` is already installed in the `jupyter-env` Conda environment you are
using for this tutorial. To try it in action, create a new notebook and name it
`Analysis.ipynb`. Add the following code to the first cell, then run it:

```python
import numpy as np
import seaborn as sns
penguins = sns.load_dataset("penguins")
```

This simply imports some python modules and loads a dataset.

Save the notebook. Now we'll add and commit the new notebook to the Git repository:

```bash
git add Analysis.ipynb
git commit -m "Add Analysis notebook"
```

So far so good. And nothing new here compared to what we've already learned
about version control. Now let's make some changes to the notebook. First we'll
replace one of the loaded modules. Update the first cell of the notebook so that
it reads:

```python
import pandas as pd
import seaborn as sns
penguins = sns.load_dataset("penguins")
```

Then create a new cell where we'll calculate the mean of each numeric value
grouped by species. In the new cell, add the following code:

```python
penguins.groupby("species").mean(numeric_only=True)
```

Run the cell and save the notebook.

Now use `git diff` to view the changes we've made to the notebook. Run:

```bash
git diff Analysis.ipynb
```

Even with very minor modifications to the notebook the diff will contain
numerous lines that are difficult to interpret. This is because the notebook not
only contains the code, but also cell metadata and output (in this case a
table produced by the second cell).

Now let's generate a more easy-to-read diff. Run:

```bash
nbdiff -s Analysis.ipynb
```

This will use the `nbdiff` tool that comes with `nbdime` to show an inline diff
of the notebook. The `-s` flag tells `nbdiff` to only show differences for the
actual code changes, ignoring changes in metadata and output. There are a number
of flags you can use here to customise the diff. The uppercase version of each
flag will ignore the respective change type. For example, to see the diff but
ignore changes to the output of cells you can run:

```bash
nbdiff -O Analysis.ipynb
```

nbdime also comes with a graphical web-based diff viewer. To try it, run:

```bash
nbdiff-web Analysis.ipynb
```

This will open up a tab in your web browser showing you changes made to the
notebook side-by-side for each cell, including also cell output. This makes it
easy to see changes made both to code and outputs such as tables and plots.

### Other tools for version control of notebooks

- You can also install the nbdime jupyter [lab
  extension](https://github.com/jupyter/nbdime) to get access to the diff
  functionality directly from the Jupyter lab interface. If you also install the
  [jupyterlab-git](https://github.com/jupyterlab/jupyterlab-git) extension you
  can both view diffs and commit changes directly from Jupyter lab.
- [VS Code](https://code.visualstudio.com/) actually comes with built-in support
  for both Jupyter notebooks and Git so that you can view [informative
  diffs](https://code.visualstudio.com/docs/datascience/jupyter-notebooks#_custom-notebook-diffing)
  directly from the editor

## Making sure notebooks work as expected

One of the great things with Jupyter notebooks is the ability to do data
exploration in an interactive way. Because loaded data, defined variables and
functions remain in the notebook until you restart the kernel, you can easily
make changes to your analysis and re-run cells to see the effect of the changes
immediately. However, this can also be a source of errors and inconsistencies if
you, during your work, modify or use variables in cells upstream of their
initial definition.

The `nbval` package can help you catch these types of errors. `nbval` is a
plugin for the `pytest` testing framework that can be used to test Jupyter
notebooks. It works by executing each cell in the notebook and comparing the
output to the output stored in the notebook. If the output is the same, the test
passes. If the output differs, the test fails. `nbval` is also pre-installed in
the `jupyter-env` Conda environment you're using for this tutorial.

As an example, we'll keep working with the `Analysis.ipynb` notebook we've
created.

Let's say we want to estimate the size of the bill of penguins using the
`bill_length_mm` and `bill_depth_mm` columns. We'll do this by adding a new cell
to our notebook with the following code:

```python
penguins["bill_size"] = (penguins["bill_length_mm"] * penguins["bill_depth_mm"])
```

Run the cell and add a new one below it. In the new cell, output the mean of
each column grouped by `island` using the following code:

```python
penguins.groupby("island").mean(numeric_only=True)
```

Run the cell to see the output. Looks good. Now we have a very simple example of
some exploratory analyses on a dataset.

Save the notebook and try running `nbval` on it to see if it works as
expected. From the commandline, run:

```bash
pytest --nbval Analysis.ipynb
```

nbval tests each cell in your notebook by executing it and comparing the output
to the output stored in the notebook. If the output is the same, the test
passes. The output of the test should look something like this:

```
collected 4 items

Analysis.ipynb ....                                                                                                   [100%]

========== 4 passed in 1.93s ==========
```

Now let's say we realize that we want to normalize the `bill_size` values by the
body mass of the penguins. We'll just modify the cell where we calculated this
value, introducing a small piece of code to divide by the `body_mass_g` column.

Change the third cell of the notebook so that it reads:

```python
penguins["bill_size"] = (penguins["bill_length_mm"] * penguins["bill_depth_mm"]) / penguins["body_mass_g"]
sns.scatterplot(data=penguins, x="bill_size", y="flipper_length_mm", hue="island")
```

Re-run the cell and save the notebook. So far so good! Let's test the notebook
again with nbval. Just like before run it from the commandline with:

```bash
pytest --nbval Analysis.ipynb
```

If you've followed the instructions, this second run of nbval should generate a
`FAILED` test, showing something like:

```
==================== short test summary info ====================
FAILED Analysis.ipynb::Cell 3
================== 1 failed, 3 passed in 1.83s ==================
```

What happened here was that we modified the cell where we calculated the
`bill_size` value, but we didn't re-run the cell where we output the mean of
each column grouped by `island`. This means that the output of the last cell in
the notebook now differs from what is actually stored in the notebook variables.
This type of error can be difficult to spot, especially if you have a large
notebook with many cells. Luckily, nbval can help us here.

> **Note** <br>
> Note that nbval reports cell numbers using 0-based numbering, so when the test
> fails on `Cell 3` it actually refers to the 4th cell in the notebook.

This problem could have been solved if we had re-run the cell where we output
the mean of each column grouped by `island`. In fact, it is good practice to
re-run all cells in a notebook before saving it. If you in addition restart the
kernel before re-running you make sure that you haven't introduced any 'hidden states'

> **Ignoring specific cells** <br>
>One caveat of `nbval` is that it doesn't work well with cells that generate
>plots. You can tell `nbval` to ignore the output of specific cells by adding
>`# NBVAL_IGNORE_OUTPUT` to the top of a cell.

> **Quick recap** <br>
> In this section we've learned:
> - How to use `nbdime` to view diffs of Jupyter notebooks
> - How to use `nbval` to test that notebooks work as expected
