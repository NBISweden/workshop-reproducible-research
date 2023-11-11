So far we've mostly worked with chunk options, which are specific to the chunk
they appear in. You can set many of these at the global document level, however,
and there are also some options specifically for tailoring the document as a
whole, regardless of chunk content.

# Document options

We've already looked at some global options, such as `title`, `author`, `format`
and `engine`. Something that would go nicely with the first two is the `date`
option. You could just write the actual date if you like, or you can use the
`today` option:

 * Add the following to the options: `date: today`

# Code folding

A useful option we haven't touched already is the `code-fold` option. This and
similar global options are specified nested inside the `format` option, like so:

```yaml
format:
    html:
        code-fold: true
```

 * Add the `code-fold` option to your document and render it.

This can be a nice default to use in scientific reports, as it hides the code by
default but is always there for those who want to inspect it. You can also use
the `code-summary` chunk option to specify a different text to show with the
folded code instead of the default `Code`, *e.g.* `code-summary: Click to show
code`.

If you want to add the `code-summary` option to all chunks you can add
the following to the yaml header:

```yaml
language:
  code-summary: Click to show code
```

You can also add the `code-tools` option, which will add a drop-down menu to
toggle visibility of all code as well as the ability to view the source of the
document.

 * Add the `code-tools: true` option and render the document.

# Table of contents

Another useful document option is to add a table of contents, which can be done
with the `toc` option. This will automatically populate the table of contents
using the headers from your document.

 * Add some more headings and/or sub-headings to your document.

 * Add the `toc: true` option to the html format and render.

The table of contents is to the right of the document by default, but you can
change it using `toc-location`. The `toc-depth` allows you to control how many
sub-heading levels are included in the table of contents.

 * Add `toc-location: left` and `toc-depth: 2` to your document and render it.

Having the table of contents on the left can be useful if you are using the
margins for something, such as we are doing in this tutorial. You can similarly
add section numbering using `number-sections` and `number-depth`. Smooth
scrolling is not enabled by default, but you can add it using `smooth-scroll:
true`. You can change the title of the table of contents using `toc-title`.

 * Add section numbers, depth, smooth scrolling and a different table of contents
  title to your document and render it.

# Themes

Quarto has a lot of [themes](https://bootswatch.com/) available for it.

 * Add `theme: flatly` under the HTML `format` option and render.

If you want to get real advanced you can play around with lots of details
regarding the themes and adjust as you see fit, or even just create your own
theme. This is a bit too advanced to go through here, but you can read about it
more in the [official documentation](https://quarto.org/docs/output-formats/html-themes.html).

# Global chunk options

The chunk options we learnt about in the previous section of this tutorial can
also be specified on the global document level. Instead of specifying *e.g.*
`warning: false` or `fig-height: 5` in individual chunks we can add it to the
main YAML header in the same manner as for *e.g.* code folding or table of
contents. We'll still have to specify options like labels or captions at the
chunk-level, though.

 * Add `warning: false` to your document header and remove it from the penguin
  figure chunk you already have.

# Embedding HTML resources

When rendering HTML documents you get any figures and other resources in a
`<document-name>_files/` directory, which is not always desirable. It's easier
to move the HTML around if all figures *etc.* are embedded directly in the HTML
itself, which can be done by specifying `embed-resources: true` in the HTML
format options. This option is false by default, meaning that you'll also have
to include the previously mentioned directory if you want to share the HTML with
anybody.

 * Remove the `<document-name>_files/` directory, refresh the rendered document
  and see what happens.

 * Add the `embed_resources` option and render your document again.

What happened first is that your figures should have disappeared when you
deleted the resources directory. Embedding resources and rendering again should
not re-create this directory, so now you'll just have a stand-alone HTML file
that is more portable than before.

# Multiple formats

So far we've mostly been working with HTML output, but you don't need to limit
yourself to a single output format if you don't want to.

 * Add the `docx: default` line in the `format:` part of your YAML header and
   render your document.

You should have gotten two separate output files now: a HTML and a DOCX (Word)
file. You can specify further options for any of the formats you include,
instead of just using the `default` settings as in this example.

 * Render your document again, but supply the `--to html` flag.

This will only render to the specified output format, which is highly useful
when you want to write a Quarto document with more than one format but not
always render them all.

# Parameters

The last document-wide option we'll touch on is *parameters*. This is useful for
when you want to be able to run the same document with different parameters or
options for some computations. How parameters are specified depends on which
engine you're using. With Knitr you can specify parameters using the `params`
option:

 * Add the following code to your YAML header:

```yaml
params:
    point_size: 2
```

 * Also change the hard-coded `geom_point(size = 2)` to `geom_point(size =
   params$point_size)` in the two `ggplot` calls in the first code chunk.

We have thus specified a parameter called `point_size` in the YAML header and
referred to it in the code using `params$point_size`. You can now change this
parameter at run-time by supplying the `-P <param>:<value>` (or `--execute-param`)
flag to `quarto render`.

Notice that this won't work if you want to use a parameter to control *e.g.* a
chunk option like `layout-ncol`. For this we need to use an in-line code
expression: `#| layout-ncol: !expr params$ncols`.

 * Add a parameter for the `layout-ncol` chunk option to the YAML header
 * Also add the `layout-ncol` chunk option to the figure chunk using the syntax
 above and render to make sure it works.

Note that to modify multiple parameters at run-time you have to use the `-P
param:value` flag multiple times, like so:

```bash
quarto render quarto-tutorial.qmd -P point_size:4 -P ncols:1
```

If you're using the Jupyter engine you can instead specify parameters by
designating a single cell as a *parameter cell*, like so:

````
```{python}
#| tags: [parameters]
point_size = 2
```
````

You can also specify parameters in a `params.yml` file and instruct quarto to use them with the `--execute-params params.yml` flag when rendering. Note that the parameters must be defined in the document (in the YAML header when using the `knitr` engine, or in a cell when using the `jupyter` engine). Pointing quarto to a `params.yml` file with `--execute-params` only overrides them when rendering.

Using parameters is extremely useful when you're using a workflow manager system (*e.g.* Snakemake or Nextflow), since you can easily specify sample-specific parameters from the command line directly from your workflow manager.

> **Quick recap** <br>
> In this sections we covered a number of document-wide options, including
> code-folding, table of contents, theming, HTML portability, using multiple
> output formats and parameters.
