--- Configuration extraction and normalization for shortcode input.
--- @module toastui._modules.config

local utils = require('./utils')

local M = {}

--- Extract document directory from current input file.
--- @return string|nil
function M.get_document_dir()
  if quarto.doc.input_file then
    return quarto.doc.input_file:match("(.+)/[^/]+$")
  end
  return nil
end

--- Read first positional shortcode argument as metadata key label.
--- @param args table
--- @return string|nil
function M.get_label(args)
  if args and #args > 0 then
    return utils.stringify(args[1])
  end
  return nil
end

--- Re-extract string metadata fields that can be lossy in generic conversion.
--- @param cfg table
--- @param raw_meta_cfg table
--- @return nil
local function hydrate_string_fields(cfg, raw_meta_cfg)
  local string_keys = { "file", "file-sep", "defaultView", "height", "timegridHeight", "date", "timeFormat" }
  for _, key in ipairs(string_keys) do
    if raw_meta_cfg[key] ~= nil then
      local s = pandoc.utils.stringify(raw_meta_cfg[key])
      if s and s ~= "" then
        cfg[key] = s
      end
    end
  end

  if raw_meta_cfg["file-sep"] ~= nil and (cfg["file-sep"] == nil or cfg["file-sep"] == "") then
    local raw = raw_meta_cfg["file-sep"]
    local t = pandoc.utils.type(raw)
    if t == "Inlines" or t == "MetaInlines" then
      local parts = {}
      for _, inline in ipairs(raw) do
        if inline.t == "Str" then
          table.insert(parts, inline.text)
        elseif inline.t == "Space" then
          table.insert(parts, " ")
        elseif inline.t == "SoftBreak" or inline.t == "LineBreak" then
          table.insert(parts, "\n")
        end
      end
      local result = table.concat(parts)
      if result ~= "" then
        cfg["file-sep"] = result
      end
    end
  end

  if raw_meta_cfg["file-sep"] ~= nil and (cfg["file-sep"] == nil or cfg["file-sep"] == "") then
    cfg["file-sep"] = "\t"
  end
end

--- Build effective shortcode config from metadata and inline kwargs.
--- @param label string|nil
--- @param meta table
--- @param kwargs table
--- @return table
function M.build_cfg(label, meta, kwargs)
  local cfg = {}
  local raw_meta_cfg = nil

  if label and meta and meta.toastui then
    local toastui_meta = meta.toastui
    if toastui_meta[label] then
      raw_meta_cfg = toastui_meta[label]
      cfg = utils.meta_to_lua(toastui_meta[label]) or {}
    end
  end

  if raw_meta_cfg then
    hydrate_string_fields(cfg, raw_meta_cfg)
  end

  for key, raw_val in pairs(kwargs) do
    local s = utils.get_kwarg(kwargs, key)
    if s ~= nil then
      cfg[key] = utils.parse_kwarg_value(s)
    end
  end

  return cfg
end

--- Build Calendar constructor options from normalized config.
--- @param cfg table
--- @return table
function M.build_options(cfg)
  local opts = {}
  local calendar_opts = {
    "defaultView", "useFormPopup", "useDetailPopup", "isReadOnly",
    "usageStatistics", "gridSelection", "timezone", "theme", "template",
    "week", "month"
  }

  for _, key in ipairs(calendar_opts) do
    if cfg[key] ~= nil then
      opts[key] = cfg[key]
    end
  end

  if opts.usageStatistics == nil then
    opts.usageStatistics = false
  end

  if opts.isReadOnly == nil then
    opts.isReadOnly = true
  end

  return opts
end

--- Read configured calendar definitions list.
--- @param cfg table
--- @return table|nil
function M.build_calendars(cfg)
  return cfg.calendars
end

--- Read configured initial date.
--- @param cfg table
--- @return string|nil
function M.initial_date(cfg)
  return cfg.date
end

--- Read and validate a detail-item selection.
--- @param cfg table
--- @param key string
--- @return table|nil,string|nil
local function detail_items(cfg, key)
  local value = cfg[key]
  if value == nil then
    return {}, nil
  end

  local items
  if type(value) == "string" then
    items = { value }
  elseif type(value) == "table" then
    local count = 0
    for index, _ in pairs(value) do
      if type(index) ~= "number" or index < 1 or index % 1 ~= 0 then
        return nil, "toastui: " .. key .. " must be a string or array"
      end
      count = count + 1
    end
    if count ~= #value then
      return nil, "toastui: " .. key .. " must be a string or array"
    end
    items = value
  else
    return nil, "toastui: " .. key .. " must be a string or array"
  end

  for _, item in ipairs(items) do
    if type(item) ~= "string" then
      return nil, "toastui: " .. key .. " must contain only strings"
    end
  end

  return items, nil
end

--- Read optional detail rows shown inside events.
--- @param cfg table
--- @return table|nil,string|nil
function M.event_detail_items(cfg)
  return detail_items(cfg, "eventDetailItems")
end

--- Read optional detail rows shown inside detail popups.
--- @param cfg table
--- @return table|nil,string|nil
function M.popup_detail_items(cfg)
  return detail_items(cfg, "popupDetailItems")
end

--- Read and validate the calendar time format.
--- @param cfg table
--- @return string|nil,string|nil
function M.time_format(cfg)
  local value = cfg.timeFormat
  if value == nil or value == "" then
    return "24h", nil
  end
  if value ~= "12h" and value ~= "24h" then
    return nil, "toastui: timeFormat must be '12h' or '24h'"
  end
  return value, nil
end

--- Determine whether navigation controls should be shown.
--- @param cfg table
--- @return boolean
function M.show_nav(cfg)
  local show_nav = cfg.navigation
  if show_nav == nil then return true end
  if show_nav == "false" or show_nav == false then
    return false
  end
  return true
end

return M
