# Tools for Reproducible Research

[![CC BY 4.0][cc-by-shield]][cc-by] [![Open in Gitpod](https://img.shields.io/badge/Gitpod-ready--to--code-FFB45B?logo=gitpod)][gitpod]

This repo contains files used in the tutorials in the NBIS / ELIXIR _Tools for
Reproducible Research_ course. The course itself lives on [GitHub
Pages](https://nbisweden.github.io/workshop-reproducible-research/), where you
can find all the relevant information.

## Students

If you are a student, you are most likely looking for the `tutorials/`
directory: this is where all the material you will need to go through the
course content is situated. You may also want to look at the non-rendered
lectures, which you can find in the `lectures/` directory.

## Contributors

The workshop website is automatically built using GitHub Actions when pushed to
the `main` branch. It is built using Quarto and Roy Francis' [Specky
template](https://github.com/royfrancis/specky/tree/main). If you want to
contribute material you should make sure it works and looks the way you want to
before you submit a pull request, which you can do by rendering the project
(with your changes) locally using the Docker image provided in this repository:

```bash
docker run -v $(pwd):/work ghcr.io/nbisweden/workshop-reproducible-research/build-website quarto render
```

The rendered website will then be rendered to the `docs/` directory, and you can
view it locally using _e.g._ `open docs/index.html` and browse your new content
from there.

## Licence

This work is licensed under a [Creative Commons Attribution 4.0 International License][cc-by].

[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg
[gitpod]: https://gitpod.io/#https://github.com/NBISweden/workshop-reproducible-research
