# Introduction

A *markup language* is a system for annotating text documents in order to
*e.g.* define formatting. HTML, if you are familiar with that, is an example of
a markup language. HTML uses tags, such as:

```html
<h1>Heading</h1>
<h2>Sub-heading</h2>
<a href="www.webpage.com">Link</a>
<ul>
  <li>List-item1</li>
  <li>List-item2</li>
  <li>List-item3</li>
</ul>
```

*Markdown* is a lightweight markup language which uses plain-text syntax in
order to be as unobtrusive as possible, so that a human can easily read it.
Some examples:

```no-highlight
# Heading

## Sub-heading

### Another deeper heading

A [link](http://example.com).

Text attributes _italic_, *italic*, **bold**, `monospace`.

Bullet list:

  * apples
  * oranges
  * pears
```

A markdown document can be converted to other formats, such as HTML or PDF, for
viewing in a browser or a PDF reader; the page you are reading right now is
written in markdown. Markdown is somewhat ill-defined, and as a consequence of
that there exists many implementations and extensions, although they share most
of the syntax. *R Markdown* is one such implementation/extension.

R Markdown documents can be used both to save and execute code (with a focus on
R) and to generate reports in various formats. This is done by mixing markdown
(as in the example above), and so-called code chunks in the same document. The
code itself, as well as the output it generates, can be included in the final
report.

R Markdown makes your analysis more reproducible by connecting your code,
figures and descriptive text. You can use it to make reproducible reports,
rather than *e.g.* copy-pasting figures into a Word document. You can also use
it as a notebook, in the same way as lab notebooks are used in a wet lab
setting (or as we us a Jupyter notebook in the [tutorial](jupyter.md)).

If you want to read more, here are some useful resources:

* A nice "Get Started" section, as a complement to this tutorial, is available
  at [RStudio.com](http://rmarkdown.rstudio.com/lesson-1.html).
* [R Markdown cheat sheet](
  https://www.rstudio.com/wp-content/uploads/2016/03/rmarkdown-cheatsheet-2.0.pdf)
  (also available from Help --> Cheatsheets in RStudio)
* [R Markdown reference guide](
  https://www.rstudio.com/wp-content/uploads/2015/03/rmarkdown-reference.pdf)
  (also available from Help --> Cheatsheets in RStudio)

This tutorial depends on files from the course GitHub repo. Take a look at the
[setup](setup.md) for instructions on how to set it up if you haven't done so
already. Place yourself in the `rmarkdown/` directory, activate your
`rmarkdown-env` Conda environment and start RStudio from the command line (type
`rstudio &`)

## Writing in R Markdown

Let's begin with starting RStudio and opening a new file (*File* --> *New File*
--> *R Markdown*). If you're using Conda you should have all the packages
needed, but install anything that RStudio prompts you to. In the window that
opens, select *Document and HTML* (which should be the default), and click *Ok*.
This will open a template R Markdown document for you. On the top is a so called
YAML header:

```yaml
---
title: "Untitled"
output:
    html_document:
        toc: true
---
```

!!! attention
    The header might look slightly different depending on your version of
    RStudio. If so, replace the default with the header above.

Here we can specify settings for the document, like the title and the output
format.

* Change the title to `My first R Markdown document`

Now, read through the rest of the template R Markdown document to get a feeling
for the format. As you can see, there are essentially three types of components
in an R Markdown document:

1. Text (written in R Markdown)
2. Code chunks (written in R or another [supported language](#r-markdown-and-other-languages))
3. The YAML header

Let's dig deeper into each of these in the following sections! But first, just
to get the flavor for things to come: press the little *Knit*-button located at
the top of the text editor panel in RStudio. This will prompt you to save the
Rmd file (do that), and generate the output file (an HTML file in this case).
It will also open up a preview of this file for you.

Some commonly used formatting written in markdown is shown below:

```no-highlight
# This is a heading

This is a paragraph.
This line-break will not appear in the output file.\
But this will (since the previous line ends with a backslash).

This is a new paragraph.

## This is a sub-heading

This is **bold text**, this is *italic text*, this is `monospaced text`,
and this is [a link](http://rmarkdown.rstudio.com/lesson-1.html).

An important feature of R Markdown is that you are allowed to use R code
inline to generate text by enclosing it with `r `. As an example: 112/67 is
equal to `r round(112/67, 2)`. You can also use multiple commands like this:
I like `r fruits <- c("apples","bananas"); paste(fruits, collapse = " and ")`!
```

The above markdown would generate something like this:

> ![](images/markdown_example.png)

Instead of reiterating information here, take a look on the first page (only
the first page!) of this [reference]( https://www.rstudio.com/wp-content/uploads/2015/03/rmarkdown-reference.pdf).
This will show you how to write more stuff in markdown and how it will be
displayed once the markdown document is converted to an output file (*e.g.*
HTML or PDF). An even more complete guide is available
[here](http://rmarkdown.rstudio.com/authoring_pandoc_markdown.html).

* Try out some of the markdown described above (and in the links) in your
  template R Markdown document! Press *Knit* to see the effect of your changes.
  Don't worry about the code chunks just yet, we'll come to that in a second.

!!! note "Quick recap"
    In this section you learned and tried out some of the basic markdown
    syntax.

## Code chunks

Enough about markdown, let's get to the fun part and include some code! Look at
the last code chunk in the template R Markdown document that you just created,
as an example:

````
```{r pressure, echo = FALSE}
plot(pressure)
```
````

The R code is surrounded by: ` ```{r}` and ` ``` `. The `r` indicates that the
code chunk contains R code (it is possible to add code chunks using other
languages, *e.g.* Python). After that comes an optional chunk name, `pressure`
in this case (this can be used to reference the code chunk as well as alleviate
debugging). Last comes chunk options, separated by commas (in this case there is
only one option: `echo = FALSE`).

!!! note "Code chunk names"
    Note that the code chunk name pressure has nothing to do with the code
    plot(pressure). In the latter case, pressure is a default R dataframe that
    is used in examples. The chunk name happened to be set to the string
    pressure as well, but could just as well have been called something else,
    *e.g.* "Plot pressure data".

Below are listed some useful chunk options related to evaluating and displaying
code chunks in the final file:

| Chunk option | Effect |
|------|-------|
| `echo = FALSE` | Prevents code, but not the results, from appearing in the finished file. This is a useful way to embed figures. |
| `include = FALSE` | Prevents both code and results from appearing in the finished file. R Markdown still runs the code in the chunk, and the results can be used by other chunks. |
| `eval = FALSE` | The code in the code chunk will not be run (but the code can be displayed in the finished file). Since the code is not evaluated, no results can be shown.|
| `results = "hide"` | Evaluate (and display) the code, but don't show the results. |
| `message = FALSE` | Prevents messages that are generated by code from appearing in the finished file. |
| `warning = FALSE` | Prevents warnings that are generated by code from appearing in the finished file. |

* Go back to your template R Markdown document in RStudio and locate the `cars`
  code chunk.
* Add the option `echo = FALSE`:

````
```{r cars, echo = FALSE}
summary(cars)
```
````

* How do you think this will affect the rendered file? Press *Knit* and check if
  you were right.
* Remove the `echo = FALSE` option and add `eval = FALSE` instead:

````
```{r cars, eval = FALSE}
summary(cars)
```
````

* How do you think this will affect the rendered file? Press *Knit* and check if
  you were right.
* Remove the `eval = FALSE` option and add `include = FALSE` instead:

````
```{r cars, include = FALSE}
summary(cars)
```
````

There are also some chunk options related to plots:

| Chunk option | Effect |
|------|-------|
| `fig.height = 9, fig.width = 6` | Set plot dimensions to 9x6 inches. (The default is 7x7.) |
| `out.height = "10cm", out.width = "8cm"` | Scale plot to 10x8 cm in the final output file. |
| `fig.cap = "This is a plot."` | Adds a figure caption.

* Go back to your template R Markdown document in RStudio and locate the
  `pressure` code chunk.
* Add the `fig.width` and `fig.height` options as below:

````
```{r pressure, echo = FALSE, fig.width = 6, fig.height = 4}
plot(pressure)
```
````

* Press *Knit* and look at the output. Can you see any differences?
* Now add a whole new code chunk to the end of the document. Give it the name
  `pressure 2` (code chunks have to have unique names, or no
  name). Add the `fig.width` and `out.width` options like this:

````
```{r pressure 2, echo = FALSE, fig.width = 9, out.width = "560px"}
plot(pressure)
```
````

* Press *Knit* and look at the output. Notice the difference between the two
  plots? In the second chunk we have first plotted a figure that is a fair bit
  larger (9 inches wide) than that in the first chunk. Next we have down-sized
  it in the final output, using the `out.width` option (where we need to use
  a size metric recognized by the output format, in this case "560px" which
  works for HTML).

Have you noticed the first chunk?

````
```{r setup, include = FALSE}
knitr::opts_chunk$set(echo = TRUE)
```
````

In this way we can set global chunk options, *i.e.* defaults for all chunks. In
this example, `echo` will always be set to `TRUE`, unless otherwise specified
in individual chunks.

!!! tip
    For more chunk options, have a look at page 2-3 of this [reference](
    https://www.rstudio.com/wp-content/uploads/2015/03/rmarkdown-reference.pdf).

It is also possible to create different types of interactive plots using
R Markdown. You can see some examples of this [here](http://www.htmlwidgets.org/showcase_networkD3.html).
If you want to try it out you can add the following code chunk to your document:

````
```{r}
library(networkD3)
data(MisLinks, MisNodes)
forceNetwork(Links = MisLinks, Nodes = MisNodes, Source = "source",
             Target = "target", Value = "value", NodeID = "name",
             Group = "group", opacity = 0.4)
```
````

!!! note "Quick recap"
    In this section you learned how to include code chunks and how to use chunk
    options to control how the output (code, results and figures) is displayed.

## The YAML header

Last but not least, we have the YAML header. Here is where you configure
general settings for the final output file, and a few other things.

The settings are written in [YAML format](https://en.wikipedia.org/wiki/YAML)
in the form *key: value*. Nested settings or sub-settings are indented with
spaces. In the template R Markdown document you can see that `html_document` is
nested under `output`, and in turn, `toc` is nested under `html_document` since
it is a setting for the HTML output. The table of contents (TOC) is
automatically compiled from the section headers (marked by #).

* Add a subsection header somewhere in your document using three `###`. Knit
  and look at how the table of contents is structured.
* Now set `toc: false` and knit again. What happened?
* A setting that works for HTML output is `toc_float: true`. Add that to your
  document (same indentation level as `toc: true`) and knit. What happened?
* In the same way, add the option `number_sections: true`. What happened?
* Do you think it looks weird with sections numbered with 0, *e.g.* 0.1? That is
  because the document does not contain any level-1-header. Add a header using
  only one `#` at the top of the document, just after the `setup` chunk. Knit
  and see what happens!

We can also set parameters in the YAML header. These are either character
strings, numerical values, or logicals, and they can be used in the R code in
the code chunks. Let's try it out:

* Add two parameters, `data` and `color`, to the YAML header. It should now
  look something like this:

```yaml
---
title: "Untitled"
output:
    html_document:
        toc: true
        toc_float: true
        number_sections: true
params:
    data: cars
    color: blue
---
```

* So now we have two parameters that we can use in the code! Modify the
  `pressure` code chunk so that it looks like this:

````
```{r pressure, fig.width = 6, fig.height = 4}
plot(get(params$data), col = params$color)
```
````

This will plot the dataset `cars` using the color `blue`, based on the
parameters we set in the YAML header.

* Knit and see what happens!

Later, we will learn how to set parameters using an external command.

We have up until now mainly been using `html_document` as an output format.
There are however a range of different available formats to choose between.
What is important to know, is that not all chunk settings work for all output
formats (this mainly regards settings related to rendering plots and figures),
and some YAML settings are specific for the given output format chosen.

* Take a look at this [gallery](http://rmarkdown.rstudio.com/gallery.html) of
  R Markdown documents to see what different kinds of output formats are
  possible to generate.
* Take a look at the last page of this [reference](
  https://www.rstudio.com/wp-content/uploads/2015/03/rmarkdown-reference.pdf)
  for a list of YAML header options, and what output formats they are available
  for.

!!! note "Quick recap"
    In this section you learned how to set document-wide settings in the YAML
    header, including document output type and user defined parameters.

## Rendering

You can render (sometimes called "knitting") reports in several different ways:

- Pressing the *Knit* button in RStudio (as we have done this far)
- Running the R command `render`: to Knit the file `my_file.Rmd` run
  `rmarkdown::render("my_file.Rmd")` in the R console.
- Running from the command line: `R -e 'rmarkdown::render("my_file.Rmd")'`

Using the `render` command, we can also set YAML header options and change
defaults (*i.e.* override those specified in the R Markdown document itself).
Here are a few useful arguments (see `?rmarkdown::render` for a full list):

- `output_format`: change output format, *e.g.* `html_document` or
  `pdf_document`
- `output_file` and `output_dir`: change directory and file name of the
  generated report file (defaults to the same name and directory as the .Rmd
  file)
- `params`: change parameter defaults. Note that only parameters already listed
  in the YAML header can be set, no new parameters can be defined

Try to use the `render` command to knit your template R Markdown document and
set the two parameters `data` and `color`. Hint: the `params` argument should
be a list, *e.g.*:

```
rmarkdown::render("my_file.Rmd", params = list(data = "cars", color = "green"))
```

### RStudio and R Markdown

You might already have noticed the various ways in which you can run code chunks
directly in RStudio:

- Place the cursor on an R command and press `CTRL + Enter` (Windows) or
  `Cmd + Enter` (Mac) to run that line in R.
- Select several R command lines and use the same keyboard shortcut as above to
  run those lines.
- To the right in each chunk there are two buttons; one for running the code in
  all chunks above the current chunk and one for running the code in the current
  chunk (depending on your layout, otherwise you can find the options in the
  *Run* drop-down).
* You can easily insert an empty chunk in your Rmd document in RStudio by
  pressing *Code* --> *Insert Chunk* in the menu.

Depending on your settings, the output of the chunk code will be displayed
inline in the Rmd document, or in RStudio's *Console* and *Plot* panels. To
customize this setting, press the cog-wheel next to the *Knit* button and select
either "Chunk Output Inline" or "Chunk Output in Console". Additionally, in the
top right in the editor panel in RStudio there is a button to toggle the
document outline. By making that visible you can click and jump between sections
(headers and named code chunks) in your R Markdown document.

## R Markdown and the case study

As you might remember from the [intro](tutorial_intro.md), we are attempting to
understand how lytic bacteriophages can be used as a future therapy for the
multiresistant bacteria MRSA (methicillin-resistant _Staphylococcus aureus_).
In this exercise we will use R Markdown to make a report in form of a
Supplementary Material HTML based on the outputs from the [Snakemake tutorial](snakemake.md).
Among the benefits of having the supplementary material (or even the full
manuscript) in R Markdown format are:

* It is fully transparent how the text, tables and figures were produced.
* If you get reviewer comments, or realize you've made a mistake somewhere, you
  can easily update your code and regenerate the document with the push of
  a button.
* By making report generation part of your workflow early on in a project, much
  of the eventual manuscript "writes itself". You no longer first have to
  finish the research part and _then_ start creating the tables and figures for
  the paper.

Before you start:

* Make sure that your working directory in R is `rmarkdown` in the course
  directory (Session > Set Working Directory).
* Open the file `rmarkdown/code/supplementary_material.Rmd`.

!!! note
    In this tutorial we have used Conda to install all the R packages we need,
    so that you get to practice how you can actually do this in projects of your
    own. You can, however, install things using `install.packages()` or
    `BiocManager::install()` as well, even though this makes it both less
    reproducible and more complicated in most cases. This is the code you will
    need to run in order to install everything without Conda:

    ```r
    BiocManager::install("ggplot2")
    BiocManager::install("reshape2")
    BiocManager::install("pheatmap")
    BiocManager::install("rtracklayer")
    BiocManager::install("GEOquery")
    install.packages("networkD3")
    ```

### Overview

* Let's start by taking a look at the YAML header at the top of the file. The
  parameters correspond to files (and sample IDs) that are generated by the
  MRSA analysis workflow (see the [Snakemake tutorial](snakemake.md)) and
  contain results that we want to include in the supplementary material
  document. We've also specified that we want to render to HTML.

```yaml
---
title: "Supplementary Materials"
output: html_document
params:
    counts_file: "results/tables/counts.tsv"
    multiqc_file: "intermediate/multiqc_general_stats.txt"
    rulegraph_file: "results/rulegraph.png"
    SRR_IDs: "SRR935090 SRR935091 SRR935092"
    GSM_IDs: "GSM1186459 GSM1186460 GSM1186461"
---
```

* From a reproducibility perspective it definitely makes sense to include
  information about who authored the document and the date it was generated.
  Add the two lines below to the YAML header. Note that we can include inline
  R code by using `` `r some_code` ``.

```
author: John Doe, Joan Dough, Jan Doh, Dyon Do
date: "`r format(Sys.time(), '%d %B, %Y')`"
```

!!! tip
    Make it a practice to keep track of all input files and add them as
    parameters rather than hard-coding them later in the R code.

Next, take a look at the `dependencies`, `read_params`, and `read_data` chunks.
They 1) load the required packages, 2) read the parameters and store them in
R objects to be used later in the code, and 3) read the data in the counts
file, the multiqc file, as well as fetch meta data from GEO. These chunks are
provided as is, and you do not need to edit them.

Below these chunks there is some markdown text that contains the Supplementary
Methods section. Note the use of section headers using `#` and `##`. Then there
is a Supplementary Tables and Figures section. This contains four code chunks,
each for a specific table or figure. Have a quick look at the code and see if
you can figure out what it does, but don't worry if you can't understand
everything.

Finally, there is a *Reproducibility* section which describes how the results in
the report can be reproduced. The `session_info` chunk prints information
regarding R version and which packages and versions that are used. We highly
encourage you to include this chunk in all your R Markdown reports: it's an
effortless way to increase reproducibility.

### Rendering options and paths

* Now that you have had a look at the R Markdown document, it is time to Knit!
  We will do this from the R terminal (rather than pressing *Knit*).

```r
rmarkdown::render("code/supplementary_material.Rmd", output_dir = "results")
```

The reason for this is that we can then redirect the output html file to be
saved in the `results/` directory.

Normally, while rendering, R code in the Rmd file will be executed using the
directory of the Rmd file as working directory (`rmarkdown/code` in this case).
However, it is good practice to write all code as if it would be executed from
the project root directory (`rmarkdown/` in this case). For instance, you can
see that we have specified the files in `params` with relative paths from the
project root directory. To set a different directory as working directory for
all chunks one modifies the knit options like this:

```r
knitr::opts_knit$set(root.dir = normalizePath('../'))
```

Here we set the working directory to the parent directory of the Rmd file
(`../`), in other words, the project root. Use this rather than `setwd()` while
working with Rmd files.

* Take a look at the output. You should find the html file in the `results`
  directory.

### Formatting tables and figures

You will probably get a good idea of the contents of the file, but the tables
look weird and the figures could be better formatted. Let's start by adjusting
the figures!

* Locate the `Setup` chunk. Here, we have already set `echo = FALSE`. Let's add
  some default figure options: `fig.height = 6, fig.width = 6, fig.align
  = 'center'`. This will make the figures slightly smaller than default and
  center them.

* Knit again, using the same R command as above. Do you notice any difference?
  Better, but still not perfect!

Let's improve the tables! We have not talked about tables before. There are
several options to print tables, here we will use the `kable` function which is
part of the `knitr` package.

* Go to the `Sample info` chunk. Replace the last line, `sample_info`, with:

```r
knitr::kable(sample_info)
```

* Knit again and look at the result. You should see a formatted table.
* The column names can be improved, and we could use a table legend. Change to
  use the following:

```r
knitr::kable(sample_info, caption = "Sample info",
             col.names = c("SRR", "GEO", "Strain", "Treatment"))
```

* Knit and check the result.
* Try to fix the table in the `QC statistics` chunk in the same manner. The
  column names are fine here so no need to change them, but add a table legend:
  "QC stats from FastQC". Knit and check your results.

Let's move on to the figures!

* Go to the `Counts barplot` chunk. To add a figure legend we have to use
  a chunk option (so not in the same way as for tables). Add the chunk option:

```r
fig.cap = "Counting statistics per sample, in terms of read counts for genes
           and reads not counted for various reasons."
```

* Knit and check the outcome!
* Next, add a figure legend to the figure in the `gene-heatmap` chunk. Here we
  can try out the possibility to add R code to generate the legend:

```r
fig.cap = paste0("Expression (log-10 counts) of genes with at least ",
                 max_cutoff, " counts in one sample and a CV>", cv_cutoff, ".")
```

This will use the `cv_cutoff` and `max_cutoff` variables to ensure that the
figure legend gives the same information as was used to generate the plot. Note
that figure legends are generated *after* the corresponding code chunk is
evaluated. This means we can use objects defined in the code chunk in the
legend.

* Knit and have a look at the results.

The heatmap still looks a bit odd. Let's play with the `fig.height` and
`out.height` options, like we did above, to scale the figure in a more
appropriate way. Add this to the chunk options: `fig.height = 10,
out.height = "22cm"`. Knit and check the results. Does it look better now?

* Now let's add a third figure! This time we will not plot a figure in R, but
  use an available image file showing the structure of the Snakemake workflow
  used to generate the inputs for this report. Add a new chunk at the end of
  the Supplementary Tables and Figures section containing this code:

```r
knitr::include_graphics(normalizePath(rulegraph_file))
```

* Also, add the chunk options:

```r
fig.cap = "A rule graph showing the different steps of the bioinformatic
           analysis that is included in the Snakemake workflow."
```

and:

```r
out.height = "11cm"
```

* Knit and check the results.

!!! note "Quick recap"
    In this section you learned some additional details for making nice
    R Markdown reports in a reproducible research project setting, including
    setting the root directory, adding tables as well as setting figure and
    table captions.

!!! note "R Markdown and Snakemake"
    It is definitely possible to render R Markdown documents as part of
    a Snakemake workflow. This is something we do for the final version of the
    MRSA project (in the Docker tutorial). In such cases it is advisable to
    manage the installation of R and required R packages through your conda
    environment file and use the `rmarkdown::render()` command from the shell
    section of your Snakemake rule (see the Snakefile in `docker/` for an
    example).

## Extra material

While the above tutorial teaches you all the basics of using R Markdown, there
is much more you can do with it, if you want to! Here we cover some extra
material if you're curious to learn more, but we don't consider this to be
a main part of the course.

### A nicer session info

While the default `sessionInfo()` command is highly useful for ending your
report with all the packages and their respective versions, it can be a bit hard
to read. There is, however, another version of the same command, which you can
find in the `devtools` package. By combining this command with the `markup`
result format you can get a more nicely formatted session information:

````
```{r Session info, echo = FALSE, results = "markup"}
devtools::session_info()
```
````

### R Markdown and Snakemake

Working with R Markdown in the context of a Snakemake workflow is something that
is highly useful for reproducibility and quite easy to get going with. An
important thing that you'll have to manage a bit more than usual is, however,
the working directory of the R Markdown document, which is something you can do
with parameters, the `root_directory`, `output_dir` and `output_file` special
variables. The following is a simple example of how you can write a Snakemake
rule for R Markdown:

```python
rule report:
    input:
        report = "report.Rmd"
    output:
        html = "results/report.html"
    params:
        outdir = "results"
    shell:
        """
        Rscript -e 'parameters <- list(root_directory = getwd(),
                                       parameter_a    = "first",
                                       parameter_b    = 10);
                    rmarkdown::render("{input.report}",
                                      params      = parameters,
                                      output_dir  = "{params.outdir}",
                                      output_file = "{output.html}")'
        """
```

### R Markdown and other languages

While R is the default and original language for any R Markdown document it
supports several others, including Python, bash, SQL, JavaScript, to name a few.
This means that you can actually get the reproducibility of R Markdown documents
when coding in languages other than just R, if your language is one of the
supported ones.

Two of the most important languages are Python and bash. While bash is supported
out-of-the-box directly and only requires you to specify `bash` instead of `r` in
the start of the code chunk, Python will additionally require you to have
installed the `reticulate` package. Not only does this allow you to code in
Python directly in your in R Markdown document, but the objects and variables you use in one
language/chunk will actually be available for the other language! You can read
more about the R Markdown Python engine [here](https://rstudio.github.io/reticulate/articles/r_markdown.html).

### R Markdown and LaTeX

This tutorial has been using HTML as the output format, as this is the most
common format that many data analyses are using. Some reasons for this include
not having to think about a page-layout (which is especially useful for
documents with many figures/plots, which is common for data analyses documents),
simplified creation (*i.e.* not having to think about special LaTeX commands for
PDF output) and fewer dependencies.

The PDF format has a lot going for it as well, however, such as being an
end-point for journal articles, books and similar documents, as well as being
much more powerful (meaning a steeper learning curve) than just HTML rendering:
you can customise almost anything you can think of. Not that HTML output is
lacking in options, it's just that LaTeX is simply more feature-rich.

Let's take an example: font sizes. This is something that is quite hard to do on
a per-chunk basis in HTML, but easy in LaTeX. You can change the font size of
all HTML chunks by using a custom CSS template, for instance, but in LaTeX you
can just set the font size to something before and after a chunk, like so:

````
\footnotesize
```{r Count to 3}
seq_len(3)
```
\normalsize
````

You could also do automatic figure captions using the LaTeX engine, meaning you
won't have to add "Figure X" to each caption manually. You can even have
separate groups of figure captions, such as one group for main article figures
and one group for supplementary figures - the same goes for tables, of course.

R Markdown uses LaTeX behind the scenes to render PDF documents, but you miss
some of the features that is inherent in LaTeX by going this route. There is,
thankfully, a different file format that you can use that more explicitly merges
R Markdown and all the functionality of LaTeX, called [Sweave](https://rpubs.com/YaRrr/SweaveIntro).

Sweave allows you to use any LaTeX command you want outside of R code chunks,
which is awesome for those of you who are already using LaTeX and want to
combine it with R Markdown. These files use the `.Rnw` extension rather than
`.Rmd`, with some additional changes, such as code chunks starting with `<<>>=`
and ending with `@`.

There is simply a plethora of tweaks, variables and parameters that you can use
for PDF rendering, but it can be quite overwhelming if you're just starting out.
We recommend using HTML for most things, especially data analysis, but PDF
certainly has its strengths - you could, for example, write a whole paper (and
its supplementary material) in R Markdown with LaTeX *ONLY*, without using
Microsoft Word or anything else! It doesn't get much more reproducible than
that.

### Presentations in R Markdown

R Markdown is not only for data analyses reports and papers, but you can also
create presentations with it! In fact, most of the presentations created for
this course were done using R Markdown.

A major difference between presentations in R Markdown and *e.g.* Microsoft
PowerPoint is the same as between any Markdown document (or LaTeX, for that
matter) and the more common Microsoft Word: the more usual Microsoft software is
"what you see is what you get", while Markdown/LaTeX doesn't show you the actual
output until you've rendered it. This difference is more pronounced when it
comes to presentations, as they are more visually heavy.

In essence, a R Markdown presentation works the same way as for a R Markdown
report, except some different formatting and output specifications. There are
a number of output formats you can use, but the one we've used for this course
(for no other reason than that we like it) is [Xaringan](https://github.com/yihui/xaringan).
You can install it from Conda (`r-xaringan`) like any other package and then 
specify the output format as `xaringan::moon_reader` in your YAML header. Slides are
separated using three dashes (`---`) while two dashes (`--`) signify slide
elements that should appear on-click.

Here is a bare-bones example of a R Markdown presentation using Xaringan:

````
---
title: "A R Markdown presentation"
output: xaringan::moon_reader
---

# R Markdown presentations

You can write text like normal, including all the normal Markdown formatting
such as *italics* or **bold**.

--

Anything can be separated into on-click appearance using double dashes.

## Sub-headers work fine as well

Just remember to separate your slides with three dashes!

---

# Use code chunks just like normal

This is especially useful for presentations of data analyses, since you don't
have to have a separate R Markdown or script to create the tables/figures and
then copy/paste them into a PowerPoint!

```{r, fig.height = 5, fig.width = 5, fig.align = "center"}
data(cars)
plot(cars)
```
````

Having said that, presentations is R Markdown can do *most* things that
PowerPoint can do, but it'll take some more effort. Getting something to look
like you want in a WYSIWYG-editor like PowerPoint is easier, since you're seeing
the output as you're making it, but it'll take more experimentation in
R Markdown. You can, however, automate a lot of things, such as by using CSS
templates that apply to each slide (including things such as font styles,
header, footers, and more) or like the above mentioned benefit of having both
code and its results already in your presentation without having to muck about
with copying and pasting figures and tables to a separate presentation.

For inspiration, we suggest you go to the `lectures/` directory of the course
Git repository. You should also have a look at the [official documentation](https://slides.yihui.org/xaringan/#1)
of Xaringan (which is itself a R Markdown-based presentation), as well as its
[several alternatives](https://rmarkdown.rstudio.com/lesson-11.html). We find
that using R Markdown for presentations does take about the same time or
slightly more compared to PowerPoint once you're used to it, but there's
a learning curve - as with everything else. Anything related to actual code and
presenting results it can be much quicker, however! A good exercise is to take
one of the presentations you have given in the past (such as for a lab meeting,
a journal club, *etc.*) and try to recreate that with R Markdown. Which method
of creating presentations you prefer is, ultimately, up to you and what the
your current end-goal is for the presentation.
