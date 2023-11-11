# Lectures

This directory contains all the course lectures, which were created using
[Quarto](https://quarto.org/). The directory is set up as a Quarto project,
so as to make rendering as simple as possible.

## Requirements

Install the required packages using:

```bash
conda env create -f environment.yml
```

## Rendering

Make sure your current directory is the `lectures/` directory and run the
following command:

```bash
quarto render
```

## Render a single lecture

To render only a single, specific lecture you can execute the following command:

```bash
quarto render <lecture>
```

For example, to render the Git lecture execute `quarto render git` while inside
the `lectures/` directory.
