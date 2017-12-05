# Introduction to Jupyter notebooks
The Jupyter Notebook is an open-source web application that allows you to create and share documents that contain code, equations, visualizations and text. The functionality is partly overlapping with R Markdown (see the [tutorial](rmarkdown)), in that they both use markdown and code chunks to generate reports that integrate results of computations with the code that generated them. Jupyter Notebook comes from the Python community while R Markdown was developed by RStudio, but you could use most common programming languages in either alternative. In practice though, it's quite common that R developers use Jupyter but probably not very common that Python developers use RStudio.

## What are Jupyter notebooks for?
An excellent question if I may say so! Some applications could be:

* Python is lacking a really good IDE for doing exploratory scientific data analysis, like RStudio or Matlab. Some people use it simply as an alternative for that.
* The community around Jupyter notebooks is large and dynamic, and there are tons of tools for sharing, displaying or interacting with notebooks.
* An early ambition with Jupyter notebooks, and its predecessor IPython notebooks, was to be analogous to the lab notebook used in a wet lab. It would allow the data scientist to document her day-to-day work and interweave results, ideas, and hypotheses with the code. From a reproducibility perspective this is one of the main advantages.
* Jupyter notebooks can be used, just as R Markdown, to provide a tighter connection between your data and your results by integrating results of computations with the code that generated them. It can also do this in an interactive way that makes it very appealing for sharing with others.

As always, the best way is to try it out yourself and decide what to use it for!

## Tell me more
* The [Jupyter project site](http://jupyter.org) contains a lot of information and inspiration.
* The [Jupyter Notebook documentation](https://jupyter-notebook.readthedocs.io/en/stable/).
* A [guide](http://ipywidgets.readthedocs.io/en/stable/index.html) to using widgets for creating interactive notebooks.

# Set up
This tutorial depends on files from the course BitBucket repo. Take a look at the [intro](index.md) for instructions on how to set it up if you haven't done so already. Then open up a terminal and go to `reproducible_research_course/jupyter`.

## Install Jupyter Notebook
If you have done the [Conda tutorial](conda.md) you should know how to define an environment and install packages using Conda. Create an environment containing the following packages from the `conda-forge` channel. Don't forget to activate the environment.

* `jupyter`: for running everything
* `nb_conda`: for integrating Conda with Jupyter Notebook
* `numpy`, `matplotlib`: for generating plots

(If you don't want to use Conda for some reason you can also install Jupyter with `pip3 install jupyter`.)

!!! note "A note on nomenclature"
    * Jupyter: a project to develop open-source software, open-standards, and services for interactive computing across dozens of programming languages. Lives at [jupyter.org](jupyter.org).
    * Jupyter Notebook: A web application that you use for creating and managing notebooks. One of the outputs of the Jupyter project.
    * Jupyter notebook: The actual `.ipynb` file that constitutes your notebook.

# Practical exercise
## The Jupyter Notebook dashboard
One thing that sets Jupyter Notebook apart from what you might be used to is that it's a web application, i.e. you edit and run your code from your browser. Ok, not quite everything, you first have to start the Jupyter Notebook server.

```no-highlight
$ jupyter notebook
[I 18:02:26.722 NotebookApp] Serving notebooks from local directory: /Users/arasmus/Documents/projects/reproducible_research_course/jupyter
[I 18:02:26.723 NotebookApp] 0 active kernels
[I 18:02:26.723 NotebookApp] The Jupyter Notebook is running at:
[I 18:02:26.723 NotebookApp] http://localhost:8888/?token=e03f10ccb40efc3c6154358593c410a139b76acf2cae785c
[I 18:02:26.723 NotebookApp] Use Control-C to stop this server and shut down all kernels (twice to skip confirmation).
[C 18:02:26.724 NotebookApp]

    Copy/paste this URL into your browser when you connect for the first time,
    to login with a token:
        http://localhost:8888/?token=e03f10ccb40efc3c6154358593c410a139b76acf2cae785c
[I 18:02:27.209 NotebookApp] Accepting one-time-token-authenticated connection from ::1
```

Jupyter Notebook probably opened up a web browser for you automatically, otherwise go to the adress specified in the message in the terminal. Note that the server is running locally (as [http://localhost:8888](http://localhost:8888)) so this does not require that you have an active internet connection. Also note that it says:

```no-highlight
Serving notebooks from local directory: /Users/arasmus/Documents/projects/reproducible_research_course/jupyter.
```

Everything you do in your Notebook session will be stored in this directory, so you don't lose any work if you shut down the server.

![alt text](jupyter_dashboard.png)

What you're looking at is the Notebook dashboard. This is where you manage your files, notebooks, and kernels. The Files tab shows the files in your directory. If you've done the other tutorials the file names should look familiar; they are the files needed for running the RNA-seq workflow in Snakemake. The Running tab keeps track of all your processes. The third tab, Clusters, is used for parallel computing and won't be discussed further in this tutorial. The Conda tab lets us control our Conda environments. Let's take a quick look at that. You can see that I'm currently in the `jupyter_exercise` environment.

![alt text](jupyter_conda.png)

Let's start by creating an empty notebook by selecting the Files tab and clicking New > Notebook > Python [conda env:jupyter_exercise]. This will open up a new tab or window looking like this:

![alt text](jupyter_empty_nb.png)

## The very basics
Jupyter notebooks are made up out of cells, and you are currently standing in the first cell in your notebook. The fact that it's green indicates that it's in "edit mode", so you can write stuff in it. Cells in Jypyter notebooks can be of two types: markdown or code.

* Markdown - These cells contain static material such as captions, text, lists, images and so on. You express this using Markdown, which is a lightweight markup language. Markdown documents can then be converted to other formats for viewing (the document you're reading now is written in Markdown and then converted to HTML). The format is discussed a little more in detail in the [R Markdown tutorial](rmarkdown.md). Jupyter Notebook uses a dialect of Markdown called GitHub Flavored Markdown, which is described [here](https://guides.github.com/features/mastering-markdown/).
* Code - These are the cells that actually do something, just as code chunks do in R Markdown. You can write code in dozens of languages and all do all kinds of clever tricks. You then run the code cell and any output the code generates, such as text or figures, will be displayed beneath the cell. We will get back to this in much more detail, but for now it's enough to understand that code cells are for executing code that is interpreted by a kernel (in this case the Python version in your Conda environment).

Let's use our first cell to create a header. Change the format from Code to Markdown in the drop-down list above the cell. Double click on the cell to enter editing mode (green frame) and input "# My notebook" ("#" is used in Markdown for header 1). Run the cell with Shift-Enter. Tada (hopefully)!

Before we continue, here are some shortcuts that can be useful. Note that they are only applicable when in command mode (blue frames). Most of them are also available from the menus.

* ++enter++: enter Edit mode
* ++escape++: enter Command mode
* ++ctrl+enter++: run the cell
* ++shift+enter++: run the cell and select the cell below
* ++alt+enter++: run the cell and insert a new cell below
* ++ctrl+s++: save the notebook
* ++tab++: for code completion or indentation
* ++m++/++y++: toggle between Markdown and Code cells
* ++d++-++d++: delete a cell
* ++a++/++b++: insert cells above/below current cell
* ++x++/++c++/++v++: cut/copy/paste cells

Now let's write some code! Since we chose a Python kernel, Python would be the native language to run in a cell. Enter `print("Hello world!")` in the second cell and run it. Note how the output is displayed below the cell. This interactive way of working is one of the things that sets Jupyter Notebook apart from RStudio and R Markdown. R Markdown is typically rendered top-to-bottom in one run, while you work *in* a Jupyter notebook in a different way. This has partly changed with newer versions of RStudio, but it's probably still how most people use the two tools. Another indication of this is that there is no (good) way to hide the code cells if you want to render your Jupyter notebook to a cleaner looking report (for a publication for example).

What **is** a Jupyter notebook? Let's look a little at the notebook we're currently working in. Jupyter Notebook saves it every minute or so, so you will already have it available. We can be a little meta and do this from within the notebook itself. We do it by running some shell commands in the third code cell instead of Python code. This very handy functionality is possible by prepending the command with `!`. Try `!ls` to list the files in the current directory.

Aha, we have a new file called `Untitled.ipynb`! Look at the first ten lines of the file by using `!head Untitled.ipynb`. Seems like it's just a plain old JSON file. Since it's a text file it's suitable for version control with for example Git. It turn out that GitHub and Jupyter notebooks are the best of friends, as we will see more of later. This switching between languages and whatever-works mentality is very prominent within the Jupyter notebook community.

Variables defined in cells become variables in the global namespace. You can therefore share information between cells. Try to define a function or variable in one cell and use it in the next. For example:

```python
def print_me(str):
  print(str)
```

and

```python
print_me("Hi!")
```

Your notebook should now look something like this.

![alt text](jupyter_basic.png)

We will not dwell on either using Markdown or Python; you can make really pretty notebooks with Markdown and you can code whatever you want with Python. Rather, we will focus on the Jupyter Notebook features that allow you to do a little more than that.

!!! note "Quick recap"
    In this section we've learnt:

    * That a Jupyter notebook consists of a series of cells, and that they can be either markdown or code cells.
    * That we execute the code in a code cell with the kernel that we chose when opening the notebook.
    * We can run shell commands by prepending them with `!`.
    * A Jupyter notebook is simply a text file in JSON format.

## Magics
Magics constitute a simple command language that significantly extends the power of Jupyter notebooks. There are two types of magics:

* Line magics - Commands are prepended by "%", and whose arguments only extend to the end of the line.
* Cell magics - Starts with `%%` and then applies to the whole cell. Must be written on the first line.

Now list all available magics with `%lsmagic` (which itself is a magic). You add a question mark to a magic to show the help (e.g. `%lsmagic?`). Some of them act as shortcuts for commonly used shell commands (`%ls`, `%cp`, `%cat`, ..). Others are useful for debugging and optimizing your code (`%timeit`, `%debug`, `%prun`, ..).

A very useful magic, in particular when using shell commands a lot in your work is `%%capture`. This will capture the stdout/stderr of any code cell and store them in a Python object. Run `%%capture?` to display the help and try to understand how it works. Try it out with either some Python code, other magics or shell commands.

??? note "Click to see one example"
    ```no-highlight
    %%capture output
    !echo "Print to stdout"
    !echo "Print to stderr" >&2
    ```

    and in another cell

    ```python
    print("stdout:" + output.stdout)
    print("stderr:" + output.stderr)
    ```


The `%%script` magic is used for specifying a program (bash, perl, ruby, ..) with which to run the code (similar to a shebang). For some languages it's possible to use these shortcuts:

* `%%ruby`
* `%%perl`
* `%%bash`
* `%%html`
* `%%latex`
* `%%R` (here you have to load the r2py extension with `%load_ext rpy2.ipython` first)

Try this out if you know any of the languages above. Otherwise you can always try to print the quadratic formula with LaTeX!
```no-highlight
\begin{array}{*{20}c} {x = \frac{{ - b \pm \sqrt {b^2 - 4ac} }}{{2a}}} & {{\rm{when}}} & {ax^2 + bx + c = 0} \\ \end{array}
```

Python's favorite library for plotting, matplotlib, has it's own magic as well: `%matplotlib`. Try out the code below, and you should hopefully get a pretty sine curve.

```python
%matplotlib inline
#Plot default curve
import numpy as np
import matplotlib.pyplot as plt
x = np.linspace(0,3*np.pi,100)
y = np.sin(x)
fig = plt.figure()
ax = fig.add_subplot(111)
line, = plt.plot(x, y, 'r-')
fig.canvas.draw()
```

!!! tip
    You can capture the output of a magic directly like this:
    ```python
    files_in_dir = %ls
    print(files_in_dir)
    ```

## Widgets and interactive plotting
Since we're typically running our notebooks in a web browser, they are quite well suited for also including more interactive elements. A typical use case would be that you want to communicate some results to a collaborator or to a wider audience, and that you would like them to be able to affect how the results are displayed. It could for example be to select which gene to plot for, or to see how some parameter value affects a clustering. Jupyter notebooks has great support for this in the form of widgets.

Widgets are eventful python objects that have a representation in the browser, often as a control like a slider, textbox, etc. Let's try to add a slider that allows us to change the frequency of the sine curve we plotted previously.

```python
%matplotlib notebook
# To use the widget framework, we need to import ipywidgets
import ipywidgets as widgets
import numpy as np
import matplotlib.pyplot as plt

# Plot default curve
x = np.linspace(0,3*np.pi,100)
y = np.sin(x)
fig = plt.figure()
ax = fig.add_subplot(111)
line, = plt.plot(x, y, 'r-')
fig.canvas.draw()

# Create and show the slider
slider = widgets.IntSlider(1, min = 0, max = 5)
display(slider)

# Define a function for modifying the line when the slider's value changes
def on_value_change(val):
    y = np.sin(x*val['new'])
    line.set_ydata(y)
    fig.canvas.draw_idle()

# Monitor for change, and send the new value to the function above on changes.
slider.observe(on_value_change, names='value')
```

!!! attention
    If you have problems getting this to run, first try with restarting the kernel (under the Kernel menu). Note that this will clear any variables you have loaded.

This is how it should look if everything works. You can set the frequency of the sine curve by moving the slider.
![alt text](jupyter_widget.png)

There are lots of widgets and the all work pretty much in the same way; you listen for some event to happen and if it does you pass the new state to some function. Here is a [list of all available widgets](http://ipywidgets.readthedocs.io/en/stable/examples/Widget%20List.html) together with documentation and examples.

IPython widgets, like we used here, is the most vanilla way of getting interactive graphs in Jupyter notebooks. Some other alternatives are:

* [plotly](https://plot.ly/python/ipython-notebook-tutorial) - actually an API to a web service that renders your graph and returns it for display in your Jupyter notebook. Generates very visually appealing graphs, but from a reproducibility perspective it's maybe not a good idea to be so reliant on a third party.
* [bokeh](https://bokeh.pydata.org/en/latest/docs/user_guide/notebook.html#userguide-notebook) - another popular tool for interactive graphs.

!!! tip
    Since you interact with Jupyter Notebook as a web server, it's very well suited for running in a Docker container (see the [tutorial](docker.md)). You could package your data, scripts and environment in a Docker image that also runs a Jupyter Notebook server. If you make this image available, say on Dockerhub, other researchers could then download it and interact with your data/code via the fancy interactive Jupyter notebooks that you have prepared for them.

## Running the MRSA workflow in a Jupyter notebook
https://nbviewer.jupyter.org/urls/bitbucket.org/scilifelab-lts/reproducible_research_course/raw/2a9c9fd304ee963915e0a63bbf49abfad07d3edd/jupyter/mrsa.ipynb
