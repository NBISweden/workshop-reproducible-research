Notebooks can be converted to various output formats such as HTML, PDF, LaTeX
*etc.* directly from the **File** -> **Download as** menu in the classic 
notebook interface, or via **File** -> **Save and Export Notebook As...** in 
the Jupyter lab interface. 

Conversion can also be performed on the command line using the `jupyter nbconvert` 
command. `nbconvert` is installed together with the `jupyter` Conda
package and is executed on the command line by running `jupyter nbconvert`. 
 
The syntax for converting a Jupyter notebook is:

```bash
jupyter nbconvert --to <FORMAT> notebook.ipynb
``` 

Here `<FORMAT>` can be any of `asciidoc`, `custom`, `html`, `latex`, `markdown`,
`notebook`, `pdf`, `python`, `rst`, `script`, `slides`. Converting to some 
output formats (*e.g.* PDF) may require you to install separate software such
as [Pandoc](https://pandoc.org/) or a **TeX** environment.

Try converting the `Untitled.ipynb` notebook that you have been working on so
far to HTML using `jupyter nbconvert`.

> **Tip** <br>
> If the plots in HTML rendered version of your notebook are not displayed 
> properly, try changing the `matplotlib_inline.backend_inline.set_matplotlib_formats('pdf', 'svg')`
> line to `matplotlib_inline.backend_inline.set_matplotlib_formats('retina')`.

`nbconvert` can also be used to run a Jupyter notebook from the command line
by running:
 
```bash
jupyter nbconvert --execute --to <FORMAT> notebook.ipynb 
```

`nbconvert` executes the cells in a notebook, captures the output and saves the
results in a new file. Try running it on the `Untitled.ipynb` notebook.

You can also specify a different output file with `--output <filename>`.

So in order to execute your `Untitled.ipynb` notebook and save it to a file 
named `report.html` you could run:

```bash
jupyter nbconvert --to html --output report.html --execute Untitled.ipynb
```

> **Quick recap** <br>
> In this section we've learned:
>
- How to convert Jupyter notebooks to various other formats
- How to use `nbconvert` to convert notebooks on the command line