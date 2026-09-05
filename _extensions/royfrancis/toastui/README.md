# TOAST UI Quarto Extension (Developer Guide)

This guide describes how the extension works internally, how files are organized,
and what TOAST UI functionality is supported.

## Purpose

The extension exposes a single shortcode, `toastui`, that renders a TOAST UI
Calendar widget in Quarto output.

Supported render targets:

- HTML (`html:js`)
- RevealJS (`revealjs`)

Other output formats are intentionally ignored.

## External Libraries

| Library | Role | Link |
|---|---|---|
| TOAST UI Calendar | Interactive calendar UI | https://github.com/nhn/tui.calendar |
| TOAST UI Calendar API docs | Option and event schema reference | https://nhn.github.io/tui.calendar/latest/ |
| Quarto | Shortcode/filter runtime and document build system | https://quarto.org/ |

## File Structure

| File | Responsibility |
|---|---|
| `toastui.lua` | Thin entrypoint. Wires modules and registers shortcode handler. |
| `_modules/dependencies.lua` | Registers JS/CSS dependencies once per document. |
| `_modules/utils.lua` | Shared helpers: metadata conversion, JSON serialization, arg parsing, path/format helpers. |
| `_modules/config.lua` | Builds effective config from YAML metadata + inline shortcode kwargs. |
| `_modules/events.lua` | Parses TSV/CSV-like files, normalizes booleans, validates required event fields. |
| `_modules/render.lua` | Produces widget HTML/JS, navigation toolbar, and calendar initialization code. |
| `toastui.css` | Toolbar and wrapper styling for generated calendars. |
| `assets/toastui-calendar.min.js` | Bundled upstream TOAST UI JavaScript. |
| `assets/toastui-calendar.min.css` | Bundled upstream TOAST UI stylesheet. |

## Runtime Flow

| Step | Module | What happens |
|---|---|---|
| 1 | `toastui.lua` | Checks output format support and exits with `pandoc.Null()` for unsupported formats. |
| 2 | `dependencies.lua` | Adds TOAST UI JS/CSS and extension CSS via Quarto HTML dependency API. |
| 3 | `config.lua` | Reads metadata by key (`toastui.<label>`), merges inline kwargs, normalizes string fields. |
| 4 | `events.lua` | Loads events from metadata or file (`file` wins), applies separator normalization, validates core fields. |
| 5 | `render.lua` | Generates container, optional nav, JS init (`new tui.Calendar(...)`), and optional `createEvents(...)`. |

## Configuration Semantics

| Source | Priority |
|---|---|
| YAML metadata (`toastui.<key>`) | Base |
| Inline shortcode kwargs | Override |

Events source precedence:

1. `events` from metadata are used if present
2. `file` overrides `events` when both are provided

## Design Choices

| Choice | Why |
|---|---|
| Modular Lua files | Makes behavior easier to test and maintain than one large script. |
| One-time dependency registration | Prevents duplicate script/style tags when multiple calendars are rendered in one page. |
| String-field hydration from raw metadata | Avoids Pandoc metadata edge cases (notably separator and inline scalar conversion). |
| File separator escape normalization (`\\t`, `\\n`) | Makes inline shortcode usage predictable and ergonomic. |
| Required event-field validation (`title`, `start`, `end`) | Catches malformed data early and emits warnings while keeping render resilient. |

## Supported TOAST UI Features

| Feature | Status | Notes |
|---|---|---|
| Calendar construction options (`defaultView`, `week`, `month`, etc.) | Supported | Passed through from extension config. |
| Calendars list (`calendars`) | Supported | Forwarded into constructor options. |
| Event data from metadata (`events`) | Supported | Accepts list of event objects. |
| Event data from text files (`file`, `file-sep`) | Supported | Header-driven parsing to objects. |
| Custom toolbar (prev/today/next + view buttons) | Supported | Controlled with `navigation`. |

## Not Supported / Out of Scope

| TOAST UI capability | Status | Notes |
|---|---|---|
| Non-HTML outputs (PDF/Docx/EPUB) | Not supported | Shortcode returns `pandoc.Null()`. |
| Advanced function callbacks in options (for example `eventFilter` function in YAML) | Limited | YAML cannot safely express arbitrary JS functions; use post-init custom JS if needed. |
| Built-in popup CSS dependencies (`tui-date-picker`, `tui-time-picker`) | Not bundled | If `useFormPopup` is enabled, upstream recommends adding picker styles yourself. |
| Live runtime API exposure (`cal` instance global access) | Not provided | Each widget is initialized in an IIFE without global handles. |

## Developer Tips

| Task | Recommendation |
|---|---|
| Add a new shortcode option | Implement in `config.lua`, then consume in `render.lua` or `events.lua`. |
| Change data parsing behavior | Edit `events.lua` (`parse_events_file` and validation flow). |
| Adjust UI controls | Update markup and JS handlers in `render.lua`, styles in `toastui.css`. |
| Debug odd metadata values | Add temporary `quarto.log.warning(...)` in `config.lua` around hydration/merge. |
| Keep compatibility | Prefer additive changes; avoid changing shortcode argument names. |
