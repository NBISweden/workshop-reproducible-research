You can render (sometimes called "knitting") reports in several different ways:

- Pressing the *Knit* button in RStudio (as we have done this far)
- Running the R command `render`: to Knit the file `my_file.Rmd` run
  `rmarkdown::render("my_file.Rmd")` in the R console.
- Running from the command line: `R -e 'rmarkdown::render("my_file.Rmd")'`

Using the `render` command, we can also set YAML header options and change
defaults (*i.e.* override those specified in the R Markdown document itself).
Here are a few useful arguments (see `?rmarkdown::render` for a full list):

- `output_format`: change output format, *e.g.* `html_document` or
  `pdf_document`
- `output_file` and `output_dir`: change directory and file name of the
  generated report file (defaults to the same name and directory as the .Rmd
  file)
- `params`: change parameter defaults. Note that only parameters already listed
  in the YAML header can be set, no new parameters can be defined

Try to use the `render` command to knit your template R Markdown document and
set the two parameters `data` and `color`. Hint: the `params` argument should
be a list, *e.g.*:

```
rmarkdown::render("my_file.Rmd", params = list(data = "cars", color = "green"))
```

You might already have noticed the various ways in which you can run code chunks
directly in RStudio:

- Place the cursor on an R command and press `CTRL + Enter` (Windows) or
  `Cmd + Enter` (Mac) to run that line in R.
- Select several R command lines and use the same keyboard shortcut as above to
  run those lines.
- To the right in each chunk there are two buttons; one for running the code in
  all chunks above the current chunk and one for running the code in the current
  chunk (depending on your layout, otherwise you can find the options in the
  *Run* drop-down).
* You can easily insert an empty chunk in your Rmd document in RStudio by
  pressing *Code* --> *Insert Chunk* in the menu.

Depending on your settings, the output of the chunk code will be displayed
inline in the Rmd document, or in RStudio's *Console* and *Plot* panels. To
customize this setting, press the cog-wheel next to the *Knit* button and select
either "Chunk Output Inline" or "Chunk Output in Console". Additionally, in the
top right in the editor panel in RStudio there is a button to toggle the
document outline. By making that visible you can click and jump between sections
(headers and named code chunks) in your R Markdown document.

!!! Success "Quick recap"
    In this section you learned how to render R Markdown documents into HTML
    documents using several different methods.
