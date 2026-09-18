# Quarto Leaflet Extension

## Purpose

The `leaflet` shortcode renders interactive Leaflet maps for HTML and revealjs outputs.

## Entry Point

- `leaflet.lua`: shortcode entrypoint that orchestrates config parsing, rendering, and dependency registration.

## Module Layout

- `_modules/metadata.lua`: metadata coercion helpers.
- `_modules/markers.lua`: marker parsing, normalization, and file loading.
- `_modules/config.lua`: config extraction from metadata/kwargs and merge logic.
- `_modules/json.lua`: JSON encoding helper for JavaScript payloads.
- `_modules/javascript.lua`: JavaScript string builder for map creation.
- `_modules/html.lua`: HTML escaping and inline error message helpers.
- `_modules/resources.lua`: one-time HTML dependency/resource registration.

## Resources

- `assets/leaflet.js`, `assets/leaflet.css`: bundled Leaflet assets.
- `assets/images/*`: Leaflet marker and layer control images.
- `quarto-leaflet.css`: extension-specific map/icon styles.

## Notes

- Extension metadata is defined in `_extension.yml`.
