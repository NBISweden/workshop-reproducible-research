The *Quarto* format (`.qmd`) is a multi-functional format, which is especially
useful for scientific coding and analyses. Quarto documents can be used both to
save and execute code as well as generating reports in various output formats.
This is done by mixing *markdown* and so-called *code chunks* in the same
document (we have [course materials for markdown](markdown) if you are
unfamiliar with this format). The code itself as well as the output it generates
can be included in the final report.

Quarto makes your analysis more reproducible by connecting your code, figures
and descriptive text. You can use it to make reproducible reports, rather than
*e.g.* copy-pasting figures into a Word document. You can also use it as a
notebook, in the same way as lab notebooks are used in a wet lab setting (or as
we utilise *Jupyter notebooks* in the tutorial after this one). Quarto itself
does not require any particular programming language to be installed - any
language you want to use can be installed separately. The currently supported
languages are R, Python, Julia and Observable. Quarto is fully compatible with
both R Markdown and Jupyter documents.

This tutorial depends on files from the course GitHub repo. Take a look at the
[setup](pre-course-setup) for instructions on how to set it up if you haven't
done so already. Place yourself in the `workshop-reproducible-research/tutorials/quarto/`
directory, activate your `quarto-env` Conda environment and start your text
editor or IDE of choice.

> **A note on R Markdown** <br>
> Quarto is an evolution of the [R Markdown](https://rmarkdown.rstudio.com/)
> format, which was previously used in this course. While R Markdown is a
> widely-used and excellent software for code and reports, Quarto is most easily
> thought of as "R Markdown 2.0". If you're familiar with R Markdown, you will
> find Quarto to be highly similar. The creators of both Quarto and R Markdown
> ([Posit](https://posit.co/)) have stated that R Markdown is not going to be
> deprecated, but most newer features will only come to Quarto. This means that
> if you've used R Markdown in the past *now* is a good time to make the switch,
> but you don't have to. You can check out the [Quarto
> website](https://quarto.org/docs/faq/rmarkdown.html) for more in-depth
> discussions regarding Quarto/R Markdown (dis-)similarities.
