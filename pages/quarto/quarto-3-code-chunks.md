Sometimes you want to add *chunk options* to the code chunks in your Quarto
documents. They are also in YAML format and are prefixed with a special type of
comment (`#|`). It can look something like this:

````
```{python}
#| echo: false
from datetime import date
print(date.today())
```
````

* Add the chunk option above to your document and render the document again.

Notice how we no longer see the code itself, just the output? This is because
the `echo` option specifies just that: whether we see the code or not. There are
a number of such chunk options that are useful to know about:

<table class="table table-hover table-condensed" border=1; style="width:600px; margin-left:auto; margin-right:auto;">
    <thead style="background-color:#DAE7F1">
        <tr>
            <td style="padding:5px; width:130px; text-align:center;"> <font size="3">
                <b> Chunk option </b>
            </td>
            <td style="padding:5px"> <font size="3">
                <b> Effect </b>
            </td>
        </tr>
    </thead>
    <tr>
        <td style="padding:5px; vertical-align:middle; text-align:center;"> <font size="3">
            `echo`
        </td>
        <td style="padding:5px"> <font size="3">
            Include the chunk code in the output.
        </td>
    </tr>
    <tr>
        <td style="padding:5px; vertical-align:middle; text-align:center;"> <font size="3">
            `eval`
        </td>
        <td style="padding:5px"> <font size="3">
            Evaluate the code chunk.
        </td>
    </tr>
    <tr>
        <td style="padding:5px; vertical-align:middle; text-align:center;"> <font size="3">
            `output`
        </td>
        <td style="padding:5px"> <font size="3">
            Include the results of executing the code in the output.
        </td>
    </tr>
    <tr>
        <td style="padding:5px; vertical-align:middle; text-align:center;"> <font size="3">
            `warning`
        </td>
        <td style="padding:5px"> <font size="3">
            Include warnings in the output.
        </td>
    </tr>
    <tr>
        <td style="padding:5px; vertical-align:middle; text-align:center;"> <font size="3">
            `error`
        </td>
        <td style="padding:5px"> <font size="3">
            Include errors in the output (note that this implies that errors
            executing code will not halt processing of the document).
        </td>
    </tr>
    <tr>
        <td style="padding:5px; vertical-align:middle; text-align:center;"> <font size="3">
            `include`
        </td>
        <td style="padding:5px"> <font size="3">
            Prevent both code and output from being included.
        </td>
    </tr>
</table>

* Check what happens if you change `echo: False` to `eval: False`.

Now the code in the code chunk is not run, which means that if you previously
added the python inline code it will no longer work because it depends on `date`
from the `datetime` module that we import in the code chunk. Remove the inline
code snippet if you added it. Then try rendering again. Now you should see the
code itself but it won't be run and therefore has no output.

# Figure options

There are also options related to figures, but for that we need to actually have
some code that produces a figure.

* Change the YAML header to use R instead of Python, remove the Python code
  chunk and replace it with the following (don't worry if you don't understand
  the R code itself, it's just as example):

````
```{r}
library("ggplot2")
library("palmerpenguins")
data(penguins, package = "palmerpenguins")
ggplot(penguins, aes(x      = bill_length_mm,
                     y      = body_mass_g,
                     colour = species)) +
    geom_point(size = 2) +
    theme_bw() +
    labs(x      = "Bill length (mm)",
         y      = "Body mass (g)",
         colour = "Species") +
    ggtitle("Penguin weight and bill length") +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_colour_manual(values = c("#c1dea0", "#85be42", "#425f21"))
```
````

When you've rendered the document you should see both the code and a figure
using the [Palmer Penguins dataset](https://allisonhorst.github.io/palmerpenguins/).
You should also see a warning along the lines of `Removed 2 rows containing
missing values`.

* Suppress the warning by adding `#| warning: false` as a chunk option and render.

There are two chunk options related to figure sizes: `fig-width` and
`fig-height` (expressed in inches). These allow you to experiment with your
figures and make them look the way you want.

* Add both the `fig-width: 10` and `fig-height: 5` chunk options and render.

> **Note** <br>
> These two chunk options are only available when using the Knitr engine, not
> for Jupyter. There is a way to set these for the whole document with Jupyter,
> though, which we'll talk more about in the next section of the tutorial.

You can also add captions and alt text using `fig-cap` and `fig-alt`,
respectively.

* Add a suitable caption and alt text to the figure and render.

If you want to place the caption in the margin of your document you can use the
`cap-location` chunk option.

* Add `cap-location: margin` to your chunk options and render.

> **Note** <br>
>
> On some quarto versions the `cap-location:` option may not work as expected.
> If you experience this, try also adding `#| label: fig-penguins` to the chunk.

# Cross-references

A convenient way to be able to refer to figures in text is by adding a figure
`label`, which will automatically add a figure number before your caption.

* Add a suitable label, *e.g.* `label: fig-penguins` to the chunk options.

Cross-references use the `@` symbol and the corresponding label. You can thus
write some markdown outside of a code chunk and refer to *e.g.* `@fig-penguins`,
as per the example here. This is extremely useful if you're writing a paper or a
report where you want to refer to figures and content in the markdown text.
Quarto even adds a clickable link to the figure itself as well!

# Sub-figures

It's also possible to create sub-figures using Quarto, instead of using whatever
plotting library that your created the figures with.

* Add the following (almost identical) code at the bottom of the chunk you
  already have:

```r
ggplot(penguins, aes(x      = bill_depth_mm,
                     y      = body_mass_g,
                     colour = species)) +
    geom_point(size = 2) +
    theme_bw() +
    labs(x      = "Bill depth (mm)",
         y      = "Body mass (g)",
         colour = "Species") +
    scale_colour_manual(values = c("#c1dea0", "#85be42", "#425f21"))
```

* Also add the following to the chunk options:

```no-highlight
#| fig-subcap:
#|     - Bill length vs. body mass
#|     - Bill depth vs. body mass
```

You should now see that we have two figures with separate sub-captions as well
as the overall figure caption we previously added. We can also control the
layout of these figures using the `layout-ncol` chunk option.

* Add a `layout-ncol: 2` chunk option and render the document.

We now have a different, two-column layout instead, but whether you prefer this
or just a one-column layout is up to you.

# Tables

Tables work much in the same way as figures. It might, in our example, be nice
to add a table with the data we previously plotted.

* Add the following code chunk to your document and render it:

````
```{r Penguin table}
#| label: tbl-penguins
#| tbl-cap: Palmer penguins bill length, width and body mass.
#| tbl-cap-location: margin
knitr::kable(
    penguins[1:10, c("species", "bill_length_mm", "bill_depth_mm", "body_mass_g")],
    col.names = c("Species", "Bill length (mm)", "Bill depth (mm)", "Body mass (g)")
)
```
````

> **Quick recap** <br>
> In this section you learned several chunk, figure and table options, how
> cross-referencing works and how to add sub-figures.
