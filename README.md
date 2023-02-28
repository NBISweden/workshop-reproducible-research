[![CC BY 4.0][cc-by-shield]][cc-by]


# Tools for Reproducible Research
---------------------------
Tools for Reproducible Research course

<img src="images/IRD.png" width="300" height="100" /> <img src="images/MIVEGEC.png" width="150" height="100" />

## Table of Contents

   * [Foreword](#foreword)
   * [Installation](#installation)
   * [Acknowledgement](#acknowledgement)
   * [License](#license)


## Foreword

The course itself lives [https://Juke34.github.io/workshop-reproducible-research/home](https://Juke34.github.io/workshop-reproducible-research/home),
where you can find all the relevant information.  
This work is based on the NBIS / ELIXIR course *Tools for
Reproducible Research* course.


## Installation

This part is for collaborators-teachers and developers

First, clone the current repository and install Mkdocs:

`pip install mkdocs`

When you are in the repository, add and/or modify your markdown tutorials in the docs directory.
The arborescence of the website menu is to setup in the `mkdocs.yml` file

For the theme:
`pip install mkdocs-material`

For the extensions:
`pip install pymdown-extensions`

For the plugins:
`pip install mkdocs-minify-plugin`
`pip installmkdocs-macros-plugin`

### Mkdocs commands for testing and building the website


* `mkdocs serve` - Start the live-reloading docs server, to test the site locally (http://127.0.0.1:8000/).
* `mkdocs gh-deploy` - Deploys the site on github pages.

* `mkdocs build` - Build the documentation site.
* `mkdocs new [dir-name]` - Create a new project.
* `mkdocs -h` - Print help message and exit.


### Project layout

    mkdocs.yml    # The configuration file.
    pages/
        home.md  # The documentation homepage.
        images   # Images used in the course materials
        ...      # Other markdown pages, images and other files.

### Welcome to MkDocs

    For full documentation visit [mkdocs.org](https://www.mkdocs.org).
    For full documentation about the [material mkdocs theme](https://squidfunk.github.io/mkdocs-material/)

##  Acknowledgement

 * NBIS - This work is based on the NBIS / ELIXIR course *Tools for Reproducible Research* course that can be find [here](https://github.com/NBISweden/workshop-reproducible-research).
 [<img align="right" src="images/NBIS.png" width="200" height="100" />](https://nbis.se)
 * @jhayer and @HadrienG about help in use and rendering with mkdocs and other features.
 * All contributors

## License

This work is licensed under a [Creative Commons Attribution 4.0 International License][cc-by].

[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg
