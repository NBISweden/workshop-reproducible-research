While the tutorial teaches you all the basics of using R Markdown, there
is much more you can do with it, if you want to! Here we cover some extra
material if you're curious to learn more, but we don't consider this to be
a main part of the course.

If you want to read more about R Markdown in general, here are some useful
resources:

* A nice "Get Started" section, as a complement to this tutorial, is available
  at [RStudio.com](http://rmarkdown.rstudio.com/lesson-1.html).
* [R Markdown cheat sheet](https://raw.githubusercontent.com/rstudio/cheatsheets/main/rmarkdown.pdf)
  (also available from Help --> Cheatsheets in RStudio)
* [R Markdown reference guide](
  https://www.rstudio.com/wp-content/uploads/2015/03/rmarkdown-reference.pdf)
  (also available from Help --> Cheatsheets in RStudio)

## A nicer session info

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

## R Markdown and workflows

Working with R Markdown in the context of a Snakemake or Nextflow workflow is something that
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

Doing it for Nextflow would look almost the same, except using Nextflow syntax
and variables.

## R Markdown and other languages

While R is the default and original language for any R Markdown document it
supports several others, including Python, bash, SQL, JavaScript, to name a few.
This means that you can actually get the reproducibility of R Markdown documents
when coding in languages other than just R, if your language is one of the
[supported ones](https://bookdown.org/yihui/rmarkdown/language-engines.html).

Two of the most important languages are Python and bash. While bash is supported
out-of-the-box directly and only requires you to specify `bash` instead of `r` in
the start of the code chunk, Python will additionally require you to have
installed the `reticulate` package. Not only does this allow you to code in
Python directly in your in R Markdown document, but the objects and variables you use in one
language/chunk will actually be available for the other language! You can read
more about the R Markdown Python engine [here](https://rstudio.github.io/reticulate/articles/r_markdown.html).

## R Markdown and LaTeX

This tutorial has been using HTML as the output format, as this is the most
common format that many data analyses are using. Some reasons for this include
not having to think about a page-layout (which is especially useful for
documents with many figures/plots, which is common for data analysis documents),
simplified creation (*i.e.* not having to think about special LaTeX commands for
PDF output) and fewer dependencies.

The PDF format has a lot going for it as well, however, such as being an
end-point for journal articles, books and similar documents, as well as being
much more powerful (meaning a steeper learning curve) than just HTML rendering:
you can customise almost anything you can think of. Not that HTML output is
lacking in options, it's just that LaTeX is more feature-rich.

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
some of the features that are inherent in LaTeX by going this route. There is,
thankfully, a different file format that you can use that more explicitly merges
R Markdown and all the functionality of LaTeX, called [Sweave](https://rpubs.com/YaRrr/SweaveIntro).

Sweave allows you to use any LaTeX command you want outside of R code chunks,
which is awesome for those of you who are already using LaTeX and want to
combine it with R Markdown. These files use the `.Rnw` extension rather than
`.Rmd`, with some additional changes, such as code chunks starting with `<<>>=`
and ending with `@`.

There is simply a plethora of tweaks, variables and parameters that you can use
for PDF rendering, but it can be quite overwhelming if you're just starting
out. We recommend using HTML for most things, especially data analysis, but PDF
certainly has its strengths - you could, for example, write a whole paper (and
its supplementary material) in R Markdown with LaTeX *only*, without using
Microsoft Word or anything else - it doesn't get much more reproducible than
that! In fact, there are several publication-ready templates on a per-journal
basis in the package [`rticles`](https://github.com/rstudio/rticles), which can
greatly facilitate this process!

## Presentations in R Markdown

R Markdown is not only useful for data analysis reports and papers, but you can also
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
specify the output format as `xaringan::moon_reader` in your YAML header.
Slides are separated using three dashes (`---`) while two dashes (`--`) signify
slide elements that should appear on-click.

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
presenting results can be much quicker, however! A good exercise is to take one
of the presentations you have given in the past (such as for a lab meeting, a
journal club, *etc.*) and try to recreate that with R Markdown. Which method of
creating presentations you prefer is, ultimately, up to you and what the your
current end-goal is for the presentation.
