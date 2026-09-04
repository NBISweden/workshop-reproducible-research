--- @module items
--- Accordion item extraction from YAML metadata or shortcode kwargs.

local function load_module(path)
  return require(quarto.utils.resolve_path(path):gsub("%.lua$", ""))
end

local utils = load_module("_modules/utils.lua")

local M = {}

--- Parse boolean-like strings from shortcode kwargs.
--- @param value string
--- @return boolean|nil
local function parse_bool(value)
  if value == "" then return nil end
  local lower = value:lower()
  if lower == "true" then return true end
  if lower == "false" then return false end
  return nil
end

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

--- Get accordion items by label from YAML metadata.
--- Supports both map form and list-of-maps form.
--- @param meta table
--- @param accordion_id string
--- @return table|nil, string|nil
function M.from_meta(meta, accordion_id)
  local meta_accordion = meta.accordion
  if meta_accordion == nil then
    return nil, string.format("'%s': No 'accordion' entry found in document yaml metadata.", accordion_id)
  end

  if type(meta_accordion[accordion_id]) == "table" then
    return meta_accordion[accordion_id], nil
  end

  for i = 1, #meta_accordion do
    local entry = meta_accordion[i]
    if type(entry) == "table" and entry[accordion_id] ~= nil then
      return entry[accordion_id], nil
    end
  end

  return nil, string.format("'%s': Accordion entry not found in yaml metadata.", accordion_id)
end

--- Get accordion items from shortcode kwargs.
--- Supports either header/body pair or an items JSON payload.
--- @param kwargs table
--- @param accordion_id string
--- @return table|nil, string|nil
function M.from_kwargs(kwargs, accordion_id)
  local header = utils.get_kwarg(kwargs, "header")
  local body = utils.get_kwarg(kwargs, "body")
  local items_json = utils.get_kwarg(kwargs, "items")

  if items_json ~= "" and (header ~= "" or body ~= "") then
    return nil, string.format("'%s': Use either 'header'/'body' or 'items', not both.", accordion_id)
  end

  if items_json ~= "" then
    if not looks_like_json(items_json) then
      return nil, string.format("'%s': Failed to parse 'items' JSON string.", accordion_id)
    end

    local ok, parsed = pcall(quarto.json.decode, items_json)
    if not ok or type(parsed) ~= "table" then
      return nil, string.format("'%s': Failed to parse 'items' JSON string.", accordion_id)
    end
    return parsed, nil
  end

  if header ~= "" or body ~= "" then
    if header == "" then
      return nil, string.format("'%s': 'header' kwarg is missing.", accordion_id)
    end
    if body == "" then
      return nil, string.format("'%s': 'body' kwarg is missing.", accordion_id)
    end

    local item = {
      header = header,
      body = body,
    }

    local collapsed = parse_bool(utils.get_kwarg(kwargs, "collapsed"))
    if collapsed ~= nil then
      item.collapsed = collapsed
    end

    local explicit_id = utils.get_kwarg(kwargs, "id")
    if explicit_id ~= "" then
      item.id = explicit_id
    end

    return { item }, nil
  end

  return nil, string.format("'%s': 'label' kwarg specified without 'header'/'body' or 'items' kwargs.", accordion_id)
end

return M