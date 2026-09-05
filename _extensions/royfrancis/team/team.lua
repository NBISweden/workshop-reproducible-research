--- @module team
--- @license MIT
--- @author Roy Francis
--- Quarto shortcode extension entrypoint for team layouts.

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

  return "fallback"
end

--- Resolve shortcode label from positional or named arguments.
--- For inline content, a deterministic label is generated from the payload.
--- @param args table
--- @param kwargs table
--- @return string|nil, string|nil
local function resolve_label(args, kwargs)
  local label = utils.get_kwarg(kwargs, "label")
  local has_label = label ~= ""
  local has_args = #args > 0
  local has_inline = items_mod.has_inline_content(kwargs)

  if has_label and has_args then
    return nil, "Use either a positional argument or the 'label' kwarg, not both."
  end

  if not has_label and not has_args and not has_inline then
    return nil, "No arguments provided. Supply a metadata label or inline team items."
  end

  local team_id
  if has_label then
    team_id = label
  elseif has_args then
    team_id = pandoc.utils.stringify(args[1])
  else
    team_id = "team-" .. utils.hash_string(items_mod.inline_signature(kwargs))
  end

  if not utils.is_valid_label(team_id) then
    return nil, string.format(
      "'%s': Label contains invalid characters. Only letters, numbers, dashes (-) and underscores (_) are allowed.",
      team_id
    )
  end

  return team_id, nil
end

--- Resolve team items from metadata or inline shortcode kwargs.
--- @param kwargs table
--- @param meta table
--- @param team_id string
--- @return table|nil, string|nil
local function resolve_items(kwargs, meta, team_id)
  if items_mod.has_inline_content(kwargs) then
    return items_mod.from_kwargs(kwargs, team_id)
  end

  return items_mod.from_meta(meta, team_id)
end

--- Render a single team shortcode instance.
--- @param args table
--- @param kwargs table
--- @param meta table
--- @return pandoc.Blocks|pandoc.RawInline|pandoc.Strong|pandoc.Null
local function render_team(args, kwargs, meta)
  local format = detect_output_mode()
  if format == "fallback" then
    return pandoc.Null()
  end

  local user_label, label_error = resolve_label(args, kwargs)
  if label_error then
    return render_mod.error_inline(label_error, format)
  end

  local team_items, items_error = resolve_items(kwargs, meta, user_label)
  if items_error then
    return render_mod.error_inline(items_error, format)
  end

  local name_class = utils.get_kwarg(kwargs, "name_class")
  local image_class = utils.get_kwarg(kwargs, "image_class")
  local description_class = utils.get_kwarg(kwargs, "description_class")

  return render_mod.html("quarto-team-" .. user_label, team_items, user_label, name_class, image_class, description_class)
end

return {
  team = render_team,
}
