# Introduction to Jupyter Notebook
The Jupyter Notebook is an open-source web application that allows you to create and share documents that contain code, equations, visualizations and text. The functionality is partly overlapping with R Markdown (see the [tutorial](rmarkdown)), in that they both use markdown and code chunks to generate reports that integrate results of computations with the code that generated them. Jupyter Notebook comes from the Python community while R Markdown was developed by RStudio, but you could use most common programming languages in either alternative. In practice though, it's quite common that R developers use Jupyter but probably not very common that Python developers use RStudio.

TODO: Some more text here. Sharing and server and git and gist and nbviever.

As always, the best way to understand how something works is to try it out.

## Tell me more
* The [Jupyter project site](http://jupyter.org) contains a lot of information and inspiration.
* The [Jupyter Notebook documentation](https://jupyter-notebook.readthedocs.io/en/stable/).

# Set up
This tutorial depends on files from the course BitBucket repo. Take a look at the [intro](index.md) for instructions on how to set it up if you haven't done so already. Then open up a terminal and go to `reproducible_research_course/jupyter`.

## Install Jupyter Notebook
If you have done the [Conda tutorial](conda.md) you should know how to define an environment and install packages using Conda. Create an environment containing the packages `jupyter` and `nb_conda` (for managing Conda environments from Jupyter) from the `conda-forge` channel. Don't forget to activate the environment.

(If you don't want to use Conda for some reason you can also install Jupyter with `pip3 install jupyter`.)

!!! note "A note on nomenclature"
    * Jupyter: a project to develop open-source software, open-standards, and services for interactive computing across dozens of programming languages. Available at [jupyter.org](jupyter.org).
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

Jupyter Notebook probably opened up a web browser for you automatically, otherwise go to the adress specified in the message in the terminal. Note that the server is running locally (as [http://localhost:8888](http://localhost:8888)) so this does not require that you have an active internet connection. Also note that it says "Serving notebooks from local directory: /Users/arasmus/Documents/projects/reproducible_research_course/jupyter". Everything you do using Notebook will be stored in this directory, so you don't lose any work if you shut down the server.

![alt text](jupyter_dashboard.png)

What you're looking at is the Notebook dashboard. This is where you manage your files, notebooks, and kernels. The Files tab shows the files in your directory. If you've done the other tutorials the file names should look familiar; they are the files needed for running the RNA-seq workflow in Snakemake. The Running tab keeps track of all your processes. The third tab, Clusters, is used for parallel computing and won't be discussed further in this tutorial. The Conda tab lets us control our Conda environments. Let's take a quick look at that. You can see that I'm currently in the `jupyter_exercise` environment.

![alt text](jupyter_conda.png)

Let's start by creating an empty notebook by selecting the Files tab and clicking New > Notebook > Python [conda env:jupyter_exercise]. This will open up a new tab or window looking like this:

![alt text](jupyter_empty_nb.png)

Jupyter notebooks are made up out of cells, and you are currently standing in the first cell in your notebook. The fact that it's green indicates that it's in "editing mode" so you can write stuff. Cells in Jypyter notebooks can be of two types: markdown or code.

* Markdown: These cells contain static material such as captions, text, lists, images and so on. You express this using Markdown, which is a lightweight markup language. Markdown documents can then be converted to other formats for viewing (the document you're reading now is written in Markdown and then converted to HTML). The format is discussed a little more in detail in the [R Markdown tutorial](rmarkdown.md). Jupyter Notebook uses a dialect of Markdown called GitHub Flavored Markdown, which is described [here](https://guides.github.com/features/mastering-markdown/).
* Code: These are the cells that actually do something, just as code chunks do in R Markdown. You can write code in dozens of languages and all do all kinds of clever tricks. You then run the code cell and any output the code generates, such as text or images, will be displayed beneath the cell. We will get back to this in much more detail, but for now it's enough to understand that code cells are for executing code that is interpreted by a kernel (in this case the Python version in your Conda environment).

Let's use our first cell to create a header. Change the format from Code to Markdown in the drop-down list above the cell. Double click on the cell to enter editing mode (green frame) and input "# My notebook" ("#" is used in Markdown for header 1). Run the cell with Shift-Enter. Tada (hopefully)!

Before we continue, here are some shortcuts that can be useful. Note that they are only applicable when in command mode (blue frames). Most of them are also available from the menus.

  * Enter key to enter Edit mode (Escape to enter Command mode)
  * Ctrl-Enter: run the cell
  * Shift-Enter: run the cell and select the cell below
  * Alt-Enter: run the cell and insert a new cell below
  * Ctrl-s: save the notebook
  * Tab key for code completion or indentation
  * m and y to toggle between Markdown and Code cells
  * d-d to delete a cell
  * z to undo deleting
  * a/b to insert cells above/below current cell
  * x/c/v to cut/copy/paste cells
  * Up/Down or k/j to select previous/next cells
  * h for help menu for keyboard shortcuts
  * Append ? for help on commands/methods, ?? to show source
