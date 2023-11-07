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
import matplotlib_inline
matplotlib_inline.backend_inline.set_matplotlib_formats("svg")
import pandas as pd
import seaborn as sns
penguins = sns.load_dataset("penguins")
sns.barplot(penguins, x="island", y="bill_depth_mm")
```

Save the notebook. Now we'll add and commit the new notebook to the Git repository:

```bash
git add Analysis.ipynb
git commit -m "Add Analysis notebook"
```

So far so good. And nothing new here compared to what we've already learned
about version control. Now let's make some changes to the notebook. We'll simply
add a `hue` argument to the plot to see differences in bill depth between male
and female penguins. Update the code cell so that the last line looks like this:

```python
sns.barplot(penguins, x="island", y="bill_depth_mm", hue="sex")
```

Then run the cell again and save the notebook.

Now we'll use `git diff` to view the changes we've made to the notebook. We
haven't yet enabled the `nbdime` functionality so this will be a default view of
the diff with `git`. Run:

```bash
git diff Analysis.ipynb
```

Even though we only made a small change to one line of code you'll see hundreds
of lines of output that are difficult to interpret. This is because the notebook
not only contains the code, but also cell metadata and output (in this case a
plot).

Now let's enable `nbdime` to see how it can help us. Run:

```bash
nbdime config-git --enable
```

This command will enable the git integration for `nbdime`. By default this
command will only enable the integration for the current repository, use the
flag `--global` to enable it for all repositories and your user account.

Now run `git diff Analysis.ipynb` again. You will still see a very long diff,
but this time you should also see some explanatory lines such as:

```
## replaced /cells/0/execution_count:
-  1
+  2

## replaced /cells/0/outputs/0/execution_count:
-  1
+  2

## inserted before /cells/0/outputs/1:
```

This was made possible by the `nbdime` integration. It understands the structure
of the notebook and can therefore provide a more informative diff.

However, to get a diff output that is even easier to read we can use the
`nbdiff` tool. Run:

```bash
nbdiff -s Analysis.ipynb
```

This will use the `nbdiff` tool that comes with `nbdime` to show an inline diff
of the notebook. In addition, the `-s` flag tells `nbdiff` to only show
differences for the actual code changes, ignoring changes in metadata and
output. There are a number of flags you can use here to customise the diff. The
uppercase version of each flag will ignore the respective change type. For
example, to see the diff but ignore changes to the output of cells you can run:

```bash
nbdiff -O Analysis.ipynb
```

nbdime also comes with a graphical web-based diff viewer. To try it, run:

```bash
nbdiff-web Analysis.ipynb
```

This will open up a tab in your web browser showing you changes made to the
notebook side-by-side for each cell, including also cell output (you may have to
click the "Trust outputs" button at the top of the page). This makes it
easy to see changes made to plots and other outputs.

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
 
As an example, create a new notebook named `exploration.ipynb` and add the following code to the first
cell to import the `seaborn` package and load the oh so useful `penguins` dataset:

```python
import seaborn as sns
df = sns.load_dataset("penguins")
```

Let's say we want to estimate the size of the bill of penguins using the
`bill_length_mm` and `bill_depth_mm` columns. We'll do this by adding a new cell
below the first one with the following code:

```python
df["bill_size"] = (df["bill_length_mm"] * df["bill_depth_mm"])
```

Run the cell and add a new one below it. In the new cell, output the mean of
each column grouped by `island` using the following code:

```python
df.groupby("island").mean(numeric_only=True)
```

Run the cell to see the output. Looks good. Now we have a very simple example of
some exploratory analyses on a dataset.

Save the notebook and try running `nbval` on it to see if it works as
expected. From the commandline, run:

```bash
pytest --nbval exploration.ipynb
```

nbval tests each cell in your notebook by executing it and comparing the output
to the output stored in the notebook. If the output is the same, the test
passes. The output of the test should look something like this:

```
collected 3 items                                                                                                              

exploration.ipynb ....                                                                                                   [100%]

========== 3 passed in 1.93s ==========
```

Now let's say we realize that we want to normalize the `bill_size` values by the
body mass of the penguins. We'll just modify the cell where we calculated this
value, introducing a small piece of code to divide by the `body_mass_g` column.

Wouldn't it also be nice to see how our estimated `bill_size` relates to the
flipper length of the penguins? Let's add a line of code to output a scatterplot
directly from the second cell where we also calculate the new value:

Change the second cell of the notebook so that it reads:

```python
df["bill_size"] = (df["bill_length_mm"] * df["bill_depth_mm"]) / df["body_mass_g"]
sns.scatterplot(data=df, x="bill_size", y="flipper_length_mm", hue="island")
```

Re-run the cell and save the notebook. So far so good! Let's test the notebook
again with nbval. Just like before run it from the commandline with:

```bash
pytest --nbval exploration.ipynb
```

If you've followed the instructions, this second run of nbval should generate a
`FAILED` test, showing something like:

```
=================================================== short test summary info ====================================================
FAILED exploration.ipynb::Cell 2
================================================= 1 failed, 2 passed in 1.83s ==================================================
```

What happened here was that we modified the cell where we calculated the 
`bill_size` value, but we didn't re-run the cell where we output the mean of
each column grouped by `island`. This means that the output of the last cell in
the notebook now differs from what is actually stored in the notebook variables.
This type of error can be difficult to spot, especially if you have a large
notebook with many cells. Luckily, nbval can help us here.

> **Note** <br>
> Note that nbval reports cell numbers using 0-based numbering, so when the test
> fails on `Cell 2` it actually refers to the third cell in the notebook.

This problem would have been solved if we had re-run the cell where we output
the mean of each column grouped by `island`. In fact, it is good practice to
re-run all cells in a notebook before saving it. If you in addition restart the
kernel before re-running you make sure that you haven't introduced any 'hidden states'
