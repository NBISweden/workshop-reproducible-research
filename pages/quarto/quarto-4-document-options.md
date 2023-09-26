So far we've mostly worked with chunk options, which are specific to the chunk
they appear in. You can set many of these at the global document level, however,
and there are also some options specifically for tailoring the document as a
whole, regardless of chunk content.

# Document options

We've already looked at some global options, such as `title`, `author`, `format`
and `engine`. Something that would go nicely with the first two is the `date`
option. You could just write the actual date if you like, or you can use inline
R code to get today's date.

* Add the following to the options: `date: "`r Sys.Date()`"`

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
the `code-summary` option to specify a different text to show with the folded
code instead of the default `Code`, *e.g.* `code-summary: Click to show code`.

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

* Add `theme: flatly` and render.

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
deleted to resources directory. Embedding resources and rendering again should
not re-create this directory, so now you'll just have a stand-alone HTML file
that is more portable than before.

> **Quick recap** <br>
> In this sections we covered a number of document-wide options, including
> code-folding, table of contents, theming and HTML portability.
