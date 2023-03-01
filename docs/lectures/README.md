# Lectures

## Requirements

Install the required packages using:

```bash
conda env create -f environment.yml
```

## Rendering in html

Lectures in **Rmarkdown** format can be rendered using the following from
the command line:

```bash
Rscript -e 'rmarkdown::render("<Rmd-file>", "xaringan::moon_reader")'
```

Lectures in **Jupyter notebook** format can be rendered using:

```bash
jupyter nbconvert <.ipynb-file> --to slides --debug --allow-chromium-download
```

## Rendering in pdf

Lectures in **Rmarkdown** format can be rendered using the following from
the command line:

```bash
Rscript -e 'library(webshot); webshot("<Rmd-file>", "<pdf-file-output>")'
```

Lectures in **Jupyter notebook** format can be rendered using:

```bash
jupyter nbconvert <.ipynb-file> --to webpdf --debug --allow-chromium-download 
```


## Render everything

To render all lectures you can use `snakemake`:

```bash
snakemake -j 1
```
