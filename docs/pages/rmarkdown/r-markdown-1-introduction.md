<iframe id="iframepdf" src="../../../lectures/rmarkdown/rmarkdown.pdf" frameborder="0" width="640" height="480" allowfullscreen="true" mozallowfullscreen="true" webkitallowfullscreen="true"></iframe>

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
that there exist many implementations and extensions, although they share most
of the syntax. *R Markdown* is one such implementation/extension.

R Markdown documents can be used both to save and execute code and to generate
reports in various formats. This is done by mixing markdown (as in the example
above), and so-called code chunks in the same document. The code itself, as
well as the output it generates, can be included in the final report.

R Markdown makes your analysis more reproducible by connecting your code,
figures and descriptive text. You can use it to make reproducible reports,
rather than *e.g.* copy-pasting figures into a Word document. You can also use
it as a notebook, in the same way as lab notebooks are used in a wet lab
setting (or as we utilise Jupyter notebooks in the tutorial after this one).

This tutorial depends on files from the course GitHub repo. Take a look at the
[setup](pre-course-setup) for instructions on how to set it up if you haven't
done so already. Place yourself in the `workshop-reproducible-research/tutorials/rmarkdown/`
directory, activate your `rmarkdown-env` Conda environment and start RStudio
from the command line (type `rstudio &`).
