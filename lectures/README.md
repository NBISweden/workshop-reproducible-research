# Lectures

## Requirements

Install the required packages using:

```bash
conda env create -f environment.yml
```

## Rendering

Lectures in **Rmarkdown** format can be rendered using the following from
the command line:

```bash
Rscript -e 'rmarkdown::render(<Rmd-file>, "xaringan::moon_reader")'
```

Lectures in **Jupyter notebook** format can be rendered using:

```bash
jupyter nbconvert --to slides --execute <.ipynb-file>
```

## Render everything

To render all lectures you can use `snakemake`:

```bash
snakemake -j 1
```