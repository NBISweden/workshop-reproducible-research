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

!!! Success "Quick recap"
    In this section you learned how to set document-wide settings in the YAML
    header, including document output type and user defined parameters.
