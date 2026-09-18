--- @module items
--- Team item extraction from YAML metadata or shortcode kwargs.

local function load_module(path)
  return require(quarto.utils.resolve_path(path):gsub("%.lua$", ""))
end

local utils = load_module("_modules/utils.lua")

local M = {}

--- Fast pre-check to avoid noisy parser errors for obviously non-JSON payloads.
--- @param text string
--- @return boolean
local function looks_like_json(text)
  local trimmed = text:match("^%s*(.-)%s*$") or ""
  local first = trimmed:sub(1, 1)
  local last = trimmed:sub(-1)

  if first == "[" and last == "]" then return true end
  if first == "{" and last == "}" then return true end
  return false
end

--- Determine whether kwargs contain inline team content.
--- @param kwargs table
--- @return boolean
function M.has_inline_content(kwargs)
  return utils.get_kwarg(kwargs, "items") ~= ""
    or utils.get_kwarg(kwargs, "name") ~= ""
    or utils.get_kwarg(kwargs, "image") ~= ""
    or utils.get_kwarg(kwargs, "description") ~= ""
    or utils.get_kwarg(kwargs, "name_url") ~= ""
    or utils.get_kwarg(kwargs, "image_url") ~= ""
end

--- Build a deterministic signature string from inline kwargs.
--- @param kwargs table
--- @return string
function M.inline_signature(kwargs)
  local parts = {
    utils.get_kwarg(kwargs, "label"),
    utils.get_kwarg(kwargs, "items"),
    utils.get_kwarg(kwargs, "name"),
    utils.get_kwarg(kwargs, "image"),
    utils.get_kwarg(kwargs, "description"),
    utils.get_kwarg(kwargs, "name_url"),
    utils.get_kwarg(kwargs, "image_url"),
    utils.get_kwarg(kwargs, "id"),
  }
  return table.concat(parts, "|")
end

--- Get team items by label from YAML metadata.
--- Supports both map form and list-of-maps form.
--- @param meta table
--- @param team_id string
--- @return table|nil, string|nil
function M.from_meta(meta, team_id)
  local meta_team = meta.team
  if meta_team == nil then
    return nil, string.format("'%s': No 'team' entry found in document yaml metadata.", team_id)
  end

  if type(meta_team[team_id]) == "table" then
    return meta_team[team_id], nil
  end

  for i = 1, #meta_team do
    local entry = meta_team[i]
    if type(entry) == "table" and entry[team_id] ~= nil then
      return entry[team_id], nil
    end
  end

  return nil, string.format("'%s': Team entry not found in yaml metadata.", team_id)
end

--- Get team items from shortcode kwargs.
--- Supports either a single inline item or an items JSON payload.
--- @param kwargs table
--- @param team_id string
--- @return table|nil, string|nil
function M.from_kwargs(kwargs, team_id)
  local items_json = utils.get_kwarg(kwargs, "items")
  local name = kwargs.name
  local image = kwargs.image
  local description = kwargs.description
  local name_url = kwargs.name_url
  local image_url = kwargs.image_url
  local explicit_id = utils.get_kwarg(kwargs, "id")

  local has_single_item = not utils.is_empty(name)
    or not utils.is_empty(image)
    or not utils.is_empty(description)
    or not utils.is_empty(name_url)
    or not utils.is_empty(image_url)
    or explicit_id ~= ""

  if items_json ~= "" and has_single_item then
    return nil, string.format("'%s': Use either single-item kwargs or the 'items' JSON payload, not both.", team_id)
  end

  if items_json ~= "" then
    if not looks_like_json(items_json) then
      return nil, string.format("'%s': Failed to parse 'items' JSON string.", team_id)
    end

    local ok, parsed = pcall(quarto.json.decode, items_json)
    if not ok or type(parsed) ~= "table" then
      return nil, string.format("'%s': Failed to parse 'items' JSON string.", team_id)
    end

    return parsed, nil
  end

  if has_single_item then
    return {
      {
        name = name,
        image = image,
        description = description,
        name_url = name_url,
        image_url = image_url,
        id = explicit_id ~= "" and explicit_id or nil,
      },
    }, nil
  end

  return nil, string.format("'%s': Provide a metadata label, inline item kwargs, or an 'items' JSON payload.", team_id)
end

return M