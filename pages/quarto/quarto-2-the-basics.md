Let's start with creating basic Quarto document that we can work with.

# Creating Quarto documents

Quarto documents are just plain text files with the `.qmd` extension. Create a
new file called *e.g.* `quarto-tutorial.qmd` and copy the following into it:

```yaml
---
title: "Untitled Quarto Document"
author: "Jane Doe"
format: html
---
```

This is a so-called *YAML header*, which is where we specify the general
settings of the document in the form of `key: value`. The `title` and `author`
are just what they sound like, while the `format` field specifies what type of
output you want the final report to be in (alternatives include `pdf`,
`revealjs` and [many others](https://quarto.org/docs/output-formats/all-formats.html)).
Here we have specified that we want HTML output, which is perhaps the most
useful for scientific computing.

* Change the title to `My first Quarto document` and the author to your name.

Let's add some actual content to the document, starting with some basic
markdown:

* Add some text into your Quarto document (including an empty line between the
  YAML header and the text), *e.g.* the following:

```
This is my first Quarto document!

# This is a header

This is where I'll soon add some *code* related to the first header.
```

Let's see what this document looks like when it's rendered into HTML by Quarto:

* Go to the command line and type `quarto render quarto-tutorial.qmd`.

> **Rendering** <br>
> If you're using *e.g.* RStudio or VSCode to edit your Quarto document you
> might have access to a *render* button, which means you don't have to run the
> above command from the command line if you prefer.

Open your new `quarto-tutorial.html` file that was created and see what it looks
like. It's only markdown content so far, so let's add some R code using a *code
chunk*:

````
```{r}
Sys.Date()
```
````

Notice that we delimit the code chunk from the rest of the document's contents
using three backticks (`\``) and specify the R language using curly brackets
(`{r}`). The code itself just prints the current date.

* Render the document again and see what it looks like.

You can also name chunks by adding it after the language:

````
```{r Show today's date}
Sys.Date()
```
````

This is useful for debugging when something has gone wrong, since it'll be
easier to see exactly which code chunk an error happened (instead of just
showing the chunk as a number).

We can also get *in-line code* using `r <R CODE>`, like so:

```
The current date is `r Sys.Date()`.
```

* Add the example above and render the document again to make sure it worked.

# Previewing documents

Quarto has a highly useful command for whem you're working on a document:
`preview`. It's essentially a live preview of the document you're working on
that will automatically render when you introduce changes to the document.

* Type `quarto preview quarto-tutorial.qmd` in the command line.

Your default web browser should now have opened a new window with your rendered
document, while your command line should say something like the following:

```no-highlight
Watching files for changes
Browse at http://localhost:4175/
```

You can't type new commands at the moment, because the Quarto Preview command is
still running - it's watching for any new changes to the Quarto document you
specified.

* Change or add some markdown text to your Quarto document, *e.g.* `This is a
  code chunk` instead of the previous text under the first header. Make sure you
  save the document.

The HTML document in your browser should have updated to reflect your newest
changes automatically. Previewing documents is great when you want to have
continuous feedback to the changes you make and can make the process of writing
more seamless, since you don't have to manually render all the time. Previewing
will still render the entire document, however, meaning that if you have some
heavy computations you might not want to re-render on every single save. For
those cases you might instead prefer to stick with manual rendering when you are
satisfied with multiple changes. You can abort a preview like any on-going
command, *e.g.* using `Ctrl-C`.

In the rest of the tutoral it's up to you whether you want to use `preview` or
not - the tutorial will just mention when it's time to render, you decide how
that's done.

# Rendering to PDF

So far we've only rendered to HTML, but sometimes you prefer a PDF. This entails
changing the `format` option in the YAML header:

* Change the format to `pdf` in the header and render your document.

You can add any raw LaTeX commands you want to your document when you're
rendering to PDF, *e.g.* `\footnotsize` to change the font size. You also have
LaTeX-specific settings, such as setting the geometry for the whole document or
specifying a citation method. While the details of LaTeX are outside the scope
of this course, it's useful to be aware of this functionality of Quarto so that
you may use it if you already know LaTeX or if you want to learn it.

# Languages

The examples so far have been using R, but we could just as easily have used
Python. All we have to do is to change our code chunk to specify `{python}` as
language and its content to be the equivalent Python code:

````
```{python}
from datetime import date
print(date.today())
```
````

* Change the code chunk to the above Python chunk instead and render your
  document again.

So far we've had Quarto automatically determine which language *engine* should
be used, which it detects through the code chunks we've written. We can also do
this explicitly by adding `engine: knitr` or `engine: jupyter` to the YAML
header.

* Explicitly add `engine: jupyter` to your YAML header and render the document.

While this is not necessary, it can be useful to explicitly set the language for
the document, as it makes it clearer from just the YAML header what language
will be used. There are also more language-related options for Quarto, but we'll
save those for later in the tutorial.

> **Quick recap** <br>
> In this section you learned how to create, edit and render basic Quarto
> documents using different languages.
