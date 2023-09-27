Quarto can also be used to create presentations in multiple formats such as
`reveal.js` (HTML), `beamer` (PDF) and `pptx` (PowerPoint) - the most powerful
of these formats by far is the first one. Creating presentations with Quarto is
quite similar to creating general Quarto documents, with some added features to
keep in mind.

# Slides

The first thing that's needed for creating a presentation is deciding what
constitutes a slide. The default is that slides are delimited by a document's
header levels.

 * Render your document using the `--to revealjs` flag and open it.

You should now have the same document we've been working on for this tutorial in
presentation format! You can step through the slides using the arrow keys, press
`F` to go into full-screen mode, `S` to view speaker notes, `M` for the menu
(you can also click in the lower left corner to get this menu) and `ESC` to go
back.

If you've followed along you should have one level-1 header (`#`) and two
level-2 headers (`##`). Notice that the level-1 header here will render as a
blank page with just the header content on it, while the level-2 headers will
render as normal slide headers. This all looks quite nice, and we didn't even
have change a thing! Disregard that the table on the last slide doesn't fit for
now, we'll get back to it later. Another method of delimiting slides is using a
horizontal rule, `---`, which allows you more fine-grained control over slides
and their content (and is especially useful if you want to have a slide without
a title).

# Divisions

There are many ways you can add presentation-specific content to your slides,
some of which you'd recognise from *e.g.* PowerPoint functionality. First, let's
fix that issue with the table that was larger than the page. The problem here is
one of *content overflow*, which can be fixed by adding a special `{.smaller}`
*div* (a division).

 * Add the `{.smaller}` div to the table header (it should read something like
   `## A table {.smaller}`) and render.

That should have automatically re-sized the table to fit into the slide. Another
way to solve this is to make slide content scrollable.

 * Change the `{.smaller}` div to a `{.scrollable}` div and render.

Instead of re-sizing the table we now get the ability to scroll down it instead;
whichever solution you prefer is up to you.

Adding divisions of various types like this is a common thing for Quarto
presentations. Another common presentation-functionality is incremental lists,
which can also be achieved with divisions. When adding a division to slide
content we specify the division's content in a manner similar to a code chunk,
like in the following example:

````
```
## Penguin species

::: {.incremental}
 - Adelie
 - Chinstrap
 - Gentoo
:::
```
````

 * Add the code above to your document and render it.

Stepping through incremental content works the same as for stepping through
slides, *i.e.* using the arrow keys.

 * Render your document to `html` instead of `revealjs`.

Notice that Quarto rendered the HTML document just fine, even though you now
have some presentation-specific code? This allows you to switch between the
formats on-demand without having much overhead or format-specific code, which is
great when you want to present your work without having to whip out a
full-fledged presentation and all the work that goes into that!

There are other useful divisions as well, including `{.notes}` (speaker notes),
`{.aside}` (additional commentary similar to footnotes), `{.footer}` (slide
footers), which you can add in the same way as we did for the incremental list
above.

 * Pick one of the above-mentioned divisions to add to your presentation and
   render it.

> **Note** <br>
> The notes and footer divisions will appear as normal Markdown text when
> rendering to HTML, while asides will appear in the margin. These divisions
> thus represent cases you might want to avoid if you want to be completely
> format-agnostic.

# Presentation options

Just like the other formats you can specify presentation-specific options at the
document-level using the YAML header. You could, for example, add the
`{.scrollable}` or `{.smaller}` div to the entire document.

 * Add the `revealjs` format to the YAML header as well as a `scrollable: true`
   option to it.

You can also specify one of the built-in themes here.

 * Add `theme: simple` to your YAML header and render.

You can find the entire list of themes at the [Quarto website](https://quarto.org/docs/presentations/revealjs/#themes).

# Multiple columns

Sometimes you'll want to have more than one column in your presentation, which
is done with the `{.columns}` and `{.column}` divisions. The former specifies
that a section with multiple columns is starting, while the second specifies
when each column starts, like so:

```no-highlight
:::: {.columns}

::: {.column}
Left column
:::

::: {.column}
Right column
:::

::::
```

 * Add multiple columns with some content to your presentation and render it.

You can also control the widths of these columns using *e.g.* `{.column width="40%"}`.

> **Note** <br>
> The `{.columns}` div also works for a normal HTML render, so it'll look the
> same regardless of whether you output as a document or a presentation.

# Fragments

We've already learnt how to get incremental lists working, but what about
general content we want to incrementally step through? This is done with the
`{.fragment}` div.

 * Add a `{.fragment}` div to some slide content and render.

Fragments are similar to "animations" from PowerPoint and come with lots of
built-in variations, *e.g.* `fade-out`, `grow`, `strike` and [several
others](https://quarto.org/docs/presentations/revealjs/advanced.html#fragment-classes).

 * Add a fragment variant to your content, *e.g.* `{.fragment .grow}` and render
   your document.

You can also control the order in which fragments appear using the
`fragment-index=<NUMBER>` option.

 * Create a new slide and add some content with a different order of appearance
   than the order of the code. If you need help or inspiration, click below.

<details>
<summary> Click to show </summary>

```
## Why Palmer Penguins?

::: {.fragment fragment-index=2}
![](https://allisonhorst.github.io/palmerpenguins/logo.png){fig-align="center"}
:::

::: {.fragment fragment-index=1}
The goal of `palmerpenguins` is to provide a good dataset for data exploration
and visualization, as an alternative to `iris.`
:::
```

</details>

> **Quick recap** <br>
> In this section we covered how to create presentations using Quarto, including
> how to add various divisions, global slide-options, multiple columns and
> fragments.
