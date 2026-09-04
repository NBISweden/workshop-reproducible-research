# Quarto Accordion Extension

## Purpose

The `accordion` shortcode renders Bootstrap accordion blocks in HTML and revealjs outputs, with non-interactive card-style fallbacks for PDF and Typst.

## Entry Point

- `accordion.lua`: shortcode entrypoint that handles output detection, argument validation, item resolution, and rendering dispatch.

## Module Layout

- `_modules/dependencies.lua`: one-time registration of HTML/revealjs/LaTeX dependencies.
- `_modules/items.lua`: extraction and validation of accordion items from YAML metadata or shortcode kwargs.
- `_modules/render.lua`: format-specific rendering and error output helpers.
- `_modules/utils.lua`: shared helpers for escaping, kwarg coercion, block conversion, and id generation.

## Notes

- Extension metadata is defined in `_extension.yml`.
