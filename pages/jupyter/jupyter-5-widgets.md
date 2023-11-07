Since we're typically running our notebooks in a web browser, they are quite
well suited for also including more interactive elements. A typical use case
could be that you want to communicate some results to a collaborator or to
a wider audience, and that you would like them to be able to modify how the
results are displayed. It could, for example, be to select which gene to plot
for, or to see how some parameter value affects a clustering. Jupyter notebooks
has great support for this in the form of *widgets*.

Widgets are eventful Python objects that have a representation in the browser,
often as a control like a slider, text box, etc. These are implemented in the
`ipywidgets` package.

The easiest way to get started with using widgets are via the `interact` and
`interactive` functions. These functions auto-generate widgets from functions
that you define, and then call those functions when you manipulate the widgets.
This might sound abstract so let's look at an example.

Let's take the scatterplot of the penguins dataset that we generated in the
previous section and add widgets that lets us choose variables to plot as well
as coloring of the points.

First we'll import the `interactive` function from `ipywidgets`, so add the
following code to a cell and run it:

```python
from ipywidgets import interactive
```

Now, in a new cell, define a function called `scatterplot` with the code to
generate the plot itself. Also add a `palette` argument to the function so that
we can specify the colour palette to use for the plot. The function should look
like this:

```python
def scatterplot(x, y, hue, palette):
    ax = sns.scatterplot(data=penguins, x=x, y=y, hue=hue, palette=palette)
```

Run the cell and create a new cell below it.

Next, we'll use the `interactive` function to generate a widget to control the
`x`, `y`, `hue` and `palette` arguments. The `interactive` function takes a
function as its first argument, and then keyword arguments for each of the
arguments in the function. The returned value is a widget which we will store in
a variable called `interactive_scatterplot`. Add the following to a cell and run
it:

```python
interactive_scatterplot = interactive(scatterplot, 
            x=["bill_length_mm","bill_depth_mm","flipper_length_mm","body_mass_g"], 
            y=["body_mass_g","bill_length_mm","bill_depth_mm","flipper_length_mm"],
            hue=["species","island","sex"],
            palette=["Set1","Set2","Dark2","Paired2"])
```

Importantly, all parameters defined in the `scatterplot` function must be given
in the `interactive` call. The `interactive_scatterplot` widget is now tied to
the `scatterplot` function. However, we still haven't displayed the widget
itself. To do that, simply add `interactive_scatterplot` to a new cell and run it:

```python
interactive_scatterplot
```

This should show the scatterplot with drop-down menus for each of the arguments.
Try changing the `x` and `y` variables to plot by selecting from the respective
drop-downs. The `hue` drop-down now lets you change what variable to use for
colouring the points and the `palette` drop-down changes the colour palette. As
you can see, the available options in the drop-downs are the ones we specified
in the `interactive` call.

Depending on the `type` of the passed argument different types of
widgets will be created by `interactive`. For instance:

- `int` or `float` arguments will generate a slider
- `bool` arguments (True/False) will generate checkbox widgets
- `list` arguments will generate a drop-down
- `str` arguments will generate a text-box

Let's add a slider to control the size of the points. In the Seaborn package
this is controlled by the `s` argument to the `scatterplot` function. Modify the
cell with your `scatterplot` function so it looks like this (remember to run the
cell in order to update the function definition):

```python
def scatterplot(x, y, hue, palette, size=50):
    ax = sns.scatterplot(data=penguins, x=x, y=y, hue=hue, palette=palette, s=size)
```

Note that we added a `size` argument to the function and supplied it to the
Seaborn scatterplot call with `s=size`. Setting `size=50` in the function
definition means that the default size of the points will be 50.

Now we need to add a slider for the `size` argument. Update the cell where we
call the `interactive` function so that it looks like this, then run it:

```python
interactive_scatterplot = interactive(scatterplot, 
            x=["bill_length_mm","bill_depth_mm","flipper_length_mm","body_mass_g"], 
            y=["bill_length_mm","bill_depth_mm","flipper_length_mm","body_mass_g"],
            hue=["species","island","sex"],
            palette=["Set1","Set2","Dark2","Paired2"],
            size=(20,100,10))
```

Here the `size` argument is defined as a
[tuple](https://docs.python.org/3/library/stdtypes.html#tuple) which sets the
minimum value of the slider to 20, the maximum value to 100 and the step size to
10. 

Finally, re-run the cell where we displayed the `interactive_scatterplot`
widget. You should now see a slider for the `size` argument (starting at 50).
Try changing the size of the points by moving the slider.

This is how it should look if everything works.

![](images/jupyter_widget.png){ width=700px }

There are lots of widgets, _e.g._:

- Drop-down menus
- Toggle buttons
- Range sliders
- File uploader

... And much, much more. Here is a [list of all available widgets](
https://ipywidgets.readthedocs.io/en/latest/examples/Widget%20List.html)
together with documentation and examples. Some of these widgets cannot be
auto-generated by `interactive`, but fear not! Instead of relying on
auto-generation we can define the widget and supply it directly to `interactive`.

To see this in practice, we'll modify the scatterplot function to display a
title and add a color picker widget that let's us set the color of the title
text.

First define the new widget by adding the following code to a cell and running
it:

```python
colorpicker = widgets.ColorPicker(
    concise=False,
    description='Title color',
    value='blue',
    disabled=False
)
```

Now change the scatterplot function so that it looks like this, then run the
cell:

```python
def scatterplot(x, y, hue, palette, size, color):
    ax = sns.scatterplot(data=penguins, x=x, y=y, hue=hue, palette=palette, s=size)
    ax.set_title("Penguin scatterplot", color=color)
```

Next, update the cell where we call the `interactive` function so that it looks
like this, then run the cell:

```python 
interactive_scatterplot = interactive(scatterplot, 
            x=["bill_length_mm","bill_depth_mm","flipper_length_mm","body_mass_g"], 
            y=["bill_length_mm","bill_depth_mm","flipper_length_mm","body_mass_g"],
            hue=["species","island","sex"],
            palette=["Set1","Set2","Dark2","Paired2"],
            size=(20, 100, 10),
            color=colorpicker)
```

Finally, re-run the cell where we displayed the `interactive_scatterplot`. The
plot should now have a title and you should see a new color picker below the
slider for the point size. Try changing the title colour by clicking on the new 
color picker. 

> **Attention!** <br>
> Note that you may have to close the colour picker once you've made your
> choice in order to make the plot update.

## Other interactive plots

Jupyter widgets, like we used here, is the most vanilla way of getting
interactive graphs in Jupyter notebooks. Some other alternatives are:

* [altair](https://altair-viz.github.io/) is a plotting library that uses
  Vega-Lite grammar which is reminiscent of ggplot2 in R. The syntax is
  different from what we've shown here, but it's very powerful once you get the
  hang of it.
* [Plotly](https://plot.ly/python/ipython-notebook-tutorial) is actually an
  API to a web service that renders your graph and returns it for display in
  your Jupyter notebook. Generates very visually appealing graphs, but from
  a reproducibility perspective it's maybe not a good idea to be so reliant on
  a third party.
* [Bokeh](https://bokeh.pydata.org/en/latest/docs/user_guide/notebook.html#userguide-notebook)
  is another popular tool for interactive graphs. Most plotting packages for
  Python are built on top of matplotlib, but Bokeh has its own library. This
  can give a steeper learning curve if you're used to the standard packages.
* [mpld3](http://mpld3.github.io) tries to integrate matplotlib with
  Javascript and the D3js package. It doesn't scale well for very large
  datasets, but it's easy to use and works quite seamlessly.

> **Quick recap** <br>
> In this section we've learned:
>
> - How to implement interactive widgets in notebooks
