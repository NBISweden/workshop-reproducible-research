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
but this time you should also see some explanatory lines such as 

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

## Making sure notebooks work as expected

- Use [nbval]() for testing notebooks and make sure they still work as expected.

## Parameterising notebooks

- Use [papermill](https://papermill.readthedocs.io/en/latest/) to parameterise
  notebooks and run them as scripts.