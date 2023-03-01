An essential feature of Jupyter Notebooks is of course the ability to visualize
data and results via plots. A full guide to plotting in Python is beyond the
scope of this course, but we'll offer a few glimpses into the plotting landscape
of Python.

First of all, Python has a library for plotting called
[matplotlib](https://matplotlib.org/stable/index.html), which comes packed with
functionality for creating high-quality plots. Below is an example of how to
generate a line plot of a sine wave.

```python
# Import packages
import numpy as np
import matplotlib.pyplot as plt
# Generate a set of evenly spaced numbers between 0 and 100
x = np.linspace(0,3*np.pi,100)
# Use the sine function to generate y-values
y = np.sin(x)
# Plot the data
line, = plt.plot(x, y, color='red', linestyle="-")
```

By default plots are rendered in the notebook as rasterized images which can
make the quality poor. To render in scalable vector graphics format use the
`set_matplotlib_formats` from the matplotlib_inline package:

```python
import matplotlib_inline
matplotlib_inline.backend_inline.set_matplotlib_formats('pdf', 'svg')
```

Now try running the code for the sine wave plot again.

## Other packages for plotting

As we mentioned Matplotlib comes with **a lot** of functionality which is great
because it allows you to create all sorts of plots and modify them exactly to
your liking. However, this can also mean that creating very basic plots might
involve a lot of cumbersome coding, when all you want is a simple bar chart!

Fortunately there are a number of Python packages that build upon matplotlib but
with a much simplified interface. One such popular package is
[seaborn](http://seaborn.pydata.org/). Below we'll see how to generate a nice
looking bar plot with error bars.

First import the seaborn package (using an abbreviated name to simplify typing):

```python
import seaborn as sns
```

Next we'll load some example data of penguins collected at the Palmer Station,
in Antarctica.

```python
penguins = sns.load_dataset("penguins")
# Look at first 5 lines of the data
penguins.head(5)
```

The most basic way to generate a bar plot of this data with seaborn is:

```python
sns.barplot(data=penguins)
```

Simple right? Yes, but maybe not very informative. Here seaborn simply
calculates the mean of all numeric variables for the penguins and plots them
with error bars representing a 95% confidence interval.

Let's say that instead we want to plot the mean value of the body mass of the
penguins at the different islands where they were examined.

```
sns.barplot(data=penguins, x="island", y="body_mass_g", ci="sd", errwidth=.5);
```

Here we specified to use values in the 'island' column as categories for the
x-axis, and values in the 'body_mass_g' column as values for the y-axis.
The barplot function of seaborn will then calculate the mean body mass for each
island and plot the bars. With `ci="sd"` we tell the function to draw the
standard deviation as error bars, instead of computing a confidence interval.
Finally `errwidth=.5` sets the linewidth of the error bars.

If we instead want to visualize the data as a scatterplot we can use the
`sns.scatterplot` function. Let's plot the body mass vs. bill length for all
penguins and color the data points by species. We'll also move the legend
outside of the plotting area and modify the x and y-axis labels:

```python
# Store the matplotlib axes containing the plot in a variable called 'ax'
ax = sns.scatterplot(data=penguins, x="bill_length_mm", y="body_mass_g",
                     hue="species")
# Modify the labels of the plot
ax.set_xlabel("Bill length (mm)")
ax.set_ylabel("Body mass (g)")
# Set legend position outside of plot
ax.legend(bbox_to_anchor=(1,1));
```

If you want to save a plot to file you can use the `plt.savefig` function. Add
the following to the bottom of the cell with the scatterplot code:

```python
plt.savefig("scatterplot.pdf", bbox_inches="tight")
```

The `bbox_inches="tight"` setting ensures that the figure is not clipped when
saved to file.

The Seaborn [website](http://seaborn.pydata.org/) contains great tutorials and
examples of other ways to plot data!

!!! Success "Quick recap"
    In this section we've learned:

    - How to generate simple plots with `matplotlib`
    - How to import and use the `seaborn` package for plotting
    - How to save plots from notebooks to a file
