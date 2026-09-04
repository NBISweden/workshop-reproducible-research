# Quarto Team Extension

## Purpose

The `team` shortcode renders team profile grids for `html` and `revealjs` outputs.

## Entry Point

- `team.lua`: shortcode entrypoint that resolves labels, chooses metadata or inline items, and dispatches rendering.

## Module Layout

- `_modules/dependencies.lua`: one-time HTML and revealjs dependency registration.
- `_modules/items.lua`: team item extraction from YAML metadata and inline shortcode kwargs.
- `_modules/render.lua`: HTML rendering and structured error output.
- `_modules/utils.lua`: shared parsing, hashing, escaping, and content conversion helpers.

## Resources

- `team.css`: base styles for HTML and revealjs outputs.
- `team-revealjs.css`: revealjs-specific layout adjustments.

## Notes

- Supports metadata-defined teams and inline configuration via `name`/`image` kwargs or `items` JSON.
- Markdown and raw HTML in `name` and `description` are parsed through Pandoc before rendering.
- Team item ids are generated deterministically unless an explicit `id` is supplied.