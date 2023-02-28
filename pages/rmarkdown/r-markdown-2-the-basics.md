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

!!! warning
    The header might look slightly different depending on your version of
    RStudio. If so, replace the default with the header above.

Here we can specify settings for the document, like the title and the output
format.

* Change the title to `My first R Markdown document`

Now, read through the rest of the template R Markdown document to get a feeling
for the format. As you can see, there are essentially three types of components
in an R Markdown document:

1. Text (written in R Markdown)
2. Code chunks (written in R or another [supported language](https://bookdown.org/yihui/rmarkdown/language-engines.html))
3. The YAML header

Let's dig deeper into each of these in the following sections! But first, just
to get the flavor for things to come: press the little *Knit*-button located at
the top of the text editor panel in RStudio. This will prompt you to save the
Rmd file (do that), and generate the output file (an HTML file in this case).
It will also open up a preview of this file for you.

Some commonly used formatting written in markdown is shown below, which you may
recognize from the [Git tutorial](git-7-working-remotely):

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

> ![](images/markdown_example.png){ width=600px }

Instead of reiterating information here, take a look on the first page (only
the first page!) of this [reference]( https://www.rstudio.com/wp-content/uploads/2015/03/rmarkdown-reference.pdf).
This will show you how to write more stuff in markdown and how it will be
displayed once the markdown document is converted to an output file (*e.g.*
HTML or PDF). An even more complete guide is available
[here](http://rmarkdown.rstudio.com/authoring_pandoc_markdown.html).

* Try out some of the markdown described above (and in the links) in your
  template R Markdown document! Press *Knit* to see the effect of your changes.
  Don't worry about the code chunks just yet, we'll come to that in a second.

!!! Success "Quick recap"
    In this section you learned and tried out some of the basic markdown
    syntax.
