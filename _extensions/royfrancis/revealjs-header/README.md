# revealjs-header extension

A Quarto extension for adding header logos to revealjs slides.

## Files

- `_extension.yml`: Extension manifest.
- `revealjs-header.lua`: Lua filter that reads metadata and injects the header structure.
- `revealjs-header.js`: Runtime script that mounts the header in the reveal container.
- `revealjs-header.css`: Header and logo styling.

## Metadata options

Set these in `format: revealjs:`.

| Key | Type | Description |
|---|---|---|
| `header-logo-left` | string | Path/URL of left logo image. |
| `header-logo-right` | string | Path/URL of right logo image. |
| `header-logo-left-url` | string | Optional clickable URL for left logo. Absolute URL or path relative to the `.qmd` file (a relative `.qmd` target, e.g. `../index.qmd`, resolves to `.html` in the rendered link). |
| `header-logo-right-url` | string | Optional clickable URL for right logo. Absolute URL or path relative to the `.qmd` file (a relative `.qmd` target, e.g. `../index.qmd`, resolves to `.html` in the rendered link). |
| `header-logo-left-height` | string | Optional CSS height for left logo (for example `2em`, `50px`) at normal window size. See scaling note below. |
| `header-logo-right-height` | string | Optional CSS height for right logo (for example `2em`, `50px`) at normal window size. See scaling note below. |

## Design notes

- The Lua filter sanitizes metadata values and ignores empty strings.
- Empty left/right containers are still rendered to keep a stable two-column layout.
- The JavaScript mount is idempotent and guarded so it does not duplicate headers.
- Lazy reveal image attributes (`data-src`) are normalized to `src` for header images.
- Logos scale with reveal.js's own zoom factor (`Reveal.getScale()`), published by the
  JS as `--qrh-scale` on `<html>` and read by the CSS `calc()` sizing rule, so they grow
  and shrink in proportion to slide content across screen sizes instead of staying a
  fixed pixel size. `header-logo-*-height` (set as `--qrh-logo-height` inline via the
  Lua filter) is the base size at `--qrh-scale: 1`, not a fixed on-screen size.

