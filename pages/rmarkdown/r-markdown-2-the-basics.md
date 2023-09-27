Let's begin with starting RStudio and opening a new file (*File* --> *New File*
--> *R Markdown*). If you're using Conda you should have all the packages
needed, but install anything that RStudio prompts you to. In the window that
opens, select *Document and HTML* (which should be the default), and click *OK*.
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

> **Attention!** <br>
> The header might look slightly different depending on your version of
> RStudio. If so, replace the default with the header above.

Here we can specify settings for the document, like the title and the output
format.

* Change the title to `My first R Markdown document`

Now, read through the rest of the template R Markdown document to get a feeling
for the format. As you can see, there are essentially three types of components
in an R Markdown document:

1. Text (written in R Markdown)
2. Code chunks (written in R or another [supported language](https://bookdown.org/yihui/rmarkdown/language-engines.html))
3. The YAML header

We will dig deeper into each of these in the following sections! But first, just
to get the flavour for things to come: press the little *Knit*-button located at
the top of the text editor panel in RStudio. This will prompt you to save the
`rmd` file (do that), and generate the output file (an HTML file in this case).
It will also open up a preview of this file for you.

In addition to the basic [markdown functionality](markdown), R Markdown also has
some useful code-specific additions:

```rmd
An important feature of R Markdown is that you are allowed to use R code
inline to generate text by enclosing it with `r `. As an example: 112/67 is
equal to `r round(112/67, 2)`. You can also use multiple commands like this:
I like `r fruits <- c("apples","bananas"); paste(fruits, collapse = " and ")`!
```

Any in-line code starting with "r " (note the space) will be evaluated as R code
- this works for the other supported languages as well Instead of reiterating
lots of information here, take a look on the first page (only the first page!)
of this [reference](https://www.rstudio.com/wp-content/uploads/2015/03/rmarkdown-reference.pdf).
This will show you how to write more stuff in R Markdown and how it will be
displayed once the document is converted to an output file (*e.g.* HTML or PDF).
An even more complete guide is available
[here](http://rmarkdown.rstudio.com/authoring_pandoc_markdown.html).

* Try out some basic markdown in your template R Markdown document!

* Press *Knit* to see the effect of your changes. Don't worry about the code
  chunks just yet, we'll come to that in a second.

> **Quick recap** <br>
> In this section you learned how to create, edit and render basic R Markdown
> documents, as well as how to write in-line code.
