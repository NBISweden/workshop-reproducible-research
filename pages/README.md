# Uploading pages to Canvas

The `upload-page-to-canvas.sh` bash script can be used to upload or update
pages on Canvas. Usage is simple, for example:

```bash
./upload-page-to-canvas.sh conda/conda-1-introduction.md
```

This will either upload or update the `conda-1-introduction.md` page to Canvas,
as applicable, using the default parameters as specified in the script. Things
you may want to change include the `COURSE_ID` parameter and the maximum width
of the resulting page.

## Page names

Canvas uses some automation for page naming by using the original filename.
A page using the `conda-2-the-basics.md` original file will be named `Conda
2 The Basics` at Canvas. You may manually edit this to change capitalisation
or add special characters and subsequent updates to the same page will keep
those changes. For example, manually editing the page in Canvas to be `Conda 2:
The basics` will keep that change after further updates, while changing it to
*e.g.* `Conda 2: Some basics` will cause a new page with the original, unaltered
name to be created.

# Formatting for Canvas

There are some key style-choices used in formatting markdown for Canvas, which
are highlighted in this section when they vary from standard markdown.

## Blockquotes as box-replacements

A blockquote is a good replacement for non-standard markdown boxes (such as
admonitions in `mkdocs`). There are some things to keep in mind, though, such as
having to use explicit breaks (`<br>`) when you want a newline for *e.g.*
headers.

Another thing is lists: you need to make sure you have a space between any text
and the start/end of lists, otherwise they will not render as proper lists.

```no-highlight
> **Note** <br>
> This is a blockquote formatted to look nice inside Canvas. It starts with
> a header and an explicit line break, followed by some text. Then comes a
> a blank line, a list and some more text.
> 
> * The first element, with some *italic text*
> * The second element, with some **bold** text
>
> Then there's more text, but not before a blank line.
```

## Hidden boxes

If you want to have hidden boxes that can be expanded on click, you can do so
with standard HTML tags: the `<details>` tag creates a hidden box with
clickable text to expand specified inside `<summary>` tags. It can be started
and ended with horizontal rulers (`***`) - it usually looks nice to end with
a ruler regardless of content, to make a clear demarcation between box content
and normal content, but only to start the box with a ruler if it starts with
text; if it starts with a code chunk you can exclude the starting ruler.

````no-highlight
<details>
<summary> Text to always show, such as "Click to expand" </summary>
***
This text will be hidden. It starts with some text rather than a code chunk, so
there's a starting ruler.


```bash
# This is a code chunk
```

There's also a ending ruler.
***
</details>
````
