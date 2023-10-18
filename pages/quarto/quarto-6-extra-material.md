The following material contains some more advanced things that you can do with
Quarto but are not really part of the core of the Quarto material. It's a mix of
various functionalities, and you don't have to go through it if you don't want
to.

If you're interested in learning more about Quarto in general, here are some
reading tips:

 - [The Quarto documentation](https://quarto.org/docs/guide/)
 - [A gallery of Quarto examples](https://quarto.org/docs/gallery/)
 - [An awesome list of Quarto content](https://github.com/mcanouil/awesome-quarto)

# Tabsets

Sometimes you'll want to present the same content in different ways, *e.g.* the
equivalent code in different languages. Look at the following toy example:

````
::: {.panel-tabset}
## R
```{.r}
words <- c("Foo", "bar")
print(paste(words), collapse = ' ')
```

## Python
```{.python}
words = ["Foo", "bar"]
print(' '.join(words))
```
:::
````

Try adding that to a document and see that you'll get a set of tabs that change
the content of the code chunk to the respective language. This is not only
useful for showing different languages, but can be used for other situations as
well. For example, you might want to run different analyses and show them in
different tabs, or even show different interactive elements in separate tabs.

# Callouts

If you're writing some sort of documentation, tutorial or just want to draw
special attention to something, *callouts* are here for you. They render as a
coloured block with a header and content. There are five types of callouts:
`note`, `tip`, `warning`, `caution`, and `important`. As with lots of Quarto
things they are specified using a division, like so:

```
::: {.callout-note}
This is a note callout.
:::
```

The different callouts come with appropriate colours by default, which you can
change in the theme. You can also have collapsible callouts by adding the
`collapse=true` option, where `true` will have the callout collapsed by default.
You can also specify titles in the same way using the `title=<TITLE>` option.

You can change the overall appearance of callouts by using the `appearance`
option or the `callout-appearance` global option. Valid values are `default`,
`simple` and `minimal`, with decreasing usage of colours and weights. You can
also suppress the callout icons using `icon=false` or `callout-icon: false` in a
similar manner.

# Mixing R and Python

Earlier in the tutorial we showed how to change the language using the `engine`
global option, but there is actually a way to use both R and Python in the same
Quarto document. This is done via the Knitr engine and the `reticulate` R
package, which allows communication between any variables and data you store in
either R or Python code chunks. While this may not be that common of a use-case,
it's still great that it's there for those that want access to it. We won't go
through the details of this works here, but you're welcome to go and check out
the [official reticulate website](https://rstudio.github.io/reticulate/) for
yourself.

If you just want to mix R and Python in a single Quarto document without the
interoperability between the languages it's a lot simpler, though. You can
either just install the `reticulate` package (`r-reticulate` in Conda/Mamba) or
add the `python.reticulate=FALSE` chunk option to the Python chunks.

# Citations

You can actually write whole articles in Quarto! For that purpose, it's also
great that you can cite things from a bibliography as well. Specifying the
bibliography file(s) is done using the `bibliography` global option; specifying
the citation style can be done using a `csl` (Citation Style Language) file and
the `csl` global option. Citation itself is similar to cross-referencing
(`@cross-ref`), but is surrounded by square brackets: `[@citation]`. You can
read more details about citations at the [Quarto
website](https://quarto.org/docs/authoring/footnotes-and-citations.html).
