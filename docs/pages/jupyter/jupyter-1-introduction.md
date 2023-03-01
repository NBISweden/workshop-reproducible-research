<iframe id="iframehtml" src="../../../lectures/jupyter/jupyter.slides.html" frameborder="0" width="640" height="480" allowfullscreen="true" mozallowfullscreen="true" webkitallowfullscreen="true"></iframe>

The Jupyter Notebook is an open-source web application that allows you to
create and share documents that contain code, equations, visualizations and
text. The functionality is partly overlapping with R Markdown (see the
[tutorial](r-markdown-1-introduction)), in that they both use markdown and code
chunks to generate reports that integrate results of computations with the code
that generated them. Jupyter Notebook comes from the Python community while
R Markdown was developed by RStudio, but you could use most common programming
languages in either alternative. In practice though, it's quite common that
R developers use Jupyter but probably not very common that Python developers
use RStudio. Some reasons to use Jupyter include:

* Python is lacking a really good IDE for doing exploratory scientific data
  analysis, like RStudio or Matlab. Some people use Jupyter simply as an
  alternative for that.
* The community around Jupyter notebooks is large and dynamic, and there are
  lots of tools for sharing, displaying or interacting with notebooks.
* An early ambition with Jupyter notebooks (and its predecessor IPython
  notebooks) was to be analogous to the lab notebook used in a wet lab. It
  would allow the data scientist to document his or her day-to-day work and
  interweave results, ideas, and hypotheses with the code. From
  a reproducibility perspective, this is one of the main advantages.
* Jupyter notebooks can be used, just like R Markdown, to provide a tighter
  connection between your data and your results by integrating results of
  computations with the code that generated them. They can also do this in an
  interactive way that makes them very appealing for sharing with others.

As always, the best way is to try it out yourself and decide what to use it
for!

This tutorial depends on files from the course GitHub repo. Take a look at the
[setup](pre-course-setup) for instructions on how to set it up if you haven't
done so already. Then open up a terminal and go to
`workshop-reproducible-research/tutorials/jupyter` and activate your
`jupyter-env` Conda environment.

!!! Info "A note on nomenclature"
    - Jupyter: a project to develop open-source software, open-standards, and
    services for interactive computing across dozens of programming
    languages. Lives at [jupyter.org](https://jupyter.org).
    - Jupyter Notebook: A web application that you use for creating and
    managing notebooks. One of the outputs of the Jupyter project.
    - Jupyter notebook: The actual `.ipynb` file that constitutes your
    notebook.
