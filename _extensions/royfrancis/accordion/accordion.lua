--- @module accordion
--- @license MIT
--- @author Roy Francis
--- Quarto shortcode extension entrypoint for accordion content.

local function load_module(path)
  return require(quarto.utils.resolve_path(path):gsub("%.lua$", ""))
end

local deps = load_module("_modules/dependencies.lua")
local items_mod = load_module("_modules/items.lua")
local render_mod = load_module("_modules/render.lua")
local utils = load_module("_modules/utils.lua")

--- Determine the current output mode and register required dependencies.
--- @return string
local function detect_output_mode()
  if quarto.doc.is_format("revealjs") then
    deps.add_revealjs_once()
    return "html"
  end

  if quarto.doc.is_format("html:js") then
    deps.add_html_once()
    return "html"
  end

  if quarto.doc.is_format("pdf") then
    deps.add_latex_once()
    return "pdf"
  end

  if quarto.doc.is_format("typst") then
    return "typst"
  end

  return "fallback"
end

--- Resolve shortcode label from positional or named arguments.
--- @param args table
--- @param kwargs table
--- @param format string
--- @return string|nil, pandoc.RawInline|pandoc.Strong|nil
local function resolve_label(args, kwargs, format)
  local label = utils.get_kwarg(kwargs, "label")
  local has_label = label ~= ""
  local has_args = #args > 0

  if has_label and has_args then
    return nil, render_mod.error_inline(
      "Use either a positional argument or named arguments (label), not both.",
      format
    )
  end

  if not has_label and not has_args then
    return nil, render_mod.error_inline(
      "No arguments provided. Provide contents either as yaml metadata (positional argument) or inline (label kwarg).",
      format
    )
  end

  local accordion_id = has_label and label or pandoc.utils.stringify(args[1])
  if not utils.is_valid_label(accordion_id) then
    return nil, render_mod.error_inline(
      string.format(
        "'%s': Label contains invalid characters. Only letters, numbers, dashes (-) and underscores (_) are allowed.",
        accordion_id
      ),
      format
    )
  end

  return accordion_id, nil
end

--- Extract items for the resolved label from kwargs or metadata.
--- @param args table
--- @param kwargs table
--- @param meta table
--- @param accordion_id string
--- @param format string
--- @return table|nil, pandoc.RawInline|pandoc.Strong|nil
local function resolve_items(args, kwargs, meta, accordion_id, format)
  local has_label_kwarg = utils.get_kwarg(kwargs, "label") ~= ""

  local accordion_items, err
  if has_label_kwarg then
    accordion_items, err = items_mod.from_kwargs(kwargs, accordion_id)
  else
    accordion_items, err = items_mod.from_meta(meta, accordion_id)
  end

  if err ~= nil then
    return nil, render_mod.error_inline(err, format)
  end

  if accordion_items == nil or #accordion_items == 0 then
    return nil, render_mod.error_inline(
      string.format("'%s': Missing 'header' and 'body'.", accordion_id),
      format
    )
  end

  return accordion_items, nil
end

--- Render a single accordion shortcode instance.
--- @param args table
--- @param kwargs table
--- @param meta table
--- @return pandoc.Blocks|pandoc.RawInline|pandoc.Strong
local function render_accordion(args, kwargs, meta)
  local format = detect_output_mode()

  local user_label, label_error = resolve_label(args, kwargs, format)
  if label_error then return label_error end

  local accordion_items, items_error = resolve_items(args, kwargs, meta, user_label, format)
  if items_error then return items_error end

  local dom_id = "quarto-accordion-" .. user_label

  if format == "html" then
    return render_mod.html(dom_id, accordion_items, user_label)
  end
  if format == "pdf" then
    return render_mod.pdf(accordion_items, user_label)
  end
  if format == "typst" then
    return render_mod.typst(accordion_items, user_label)
  end
  return render_mod.fallback(accordion_items, user_label)
end

return {
  accordion = render_accordion,
}
