--- @module config
--- Leaflet config extraction and merge logic.

local function load_module(path)
  return require(quarto.utils.resolve_path(path):gsub("%.lua$", ""))
end

local meta = load_module("_modules/metadata.lua")
local markers = load_module("_modules/markers.lua")

local M = {}

--- Keys treated specially and not passed through to L.map options.
M.SPECIAL = {
  center = true,
  zoom = true,
  height = true,
  width = true,
  markers = true,
  ["markers-file"] = true,
  markers_file = true,
  ["markers-sep"] = true,
  markers_sep = true,
  tile = true,
}

--- Determine whether a metadata value behaves as a map-like object.
--- @param val any
--- @return boolean
local function is_meta_map_value(val)
  if type(val) ~= "table" then return false end
  if val.t ~= nil then return val.t == "MetaMap" end
  for k, _ in pairs(val) do
    if type(k) ~= "number" then return true end
  end
  return false
end

--- Build extension-level default options from leaflet metadata.
--- @param leaflet_meta table|nil
--- @return table
function M.build_leaflet_defaults_meta(leaflet_meta)
  if type(leaflet_meta) ~= "table" then return {} end

  local defaults = {}
  for k, v in pairs(leaflet_meta) do
    if k ~= "markers" then
      if M.SPECIAL[k] then
        defaults[k] = v
      elseif not is_meta_map_value(v) then
        defaults[k] = v
      end
    end
  end

  return defaults
end

--- Merge two configs while preserving marker accumulation and tile options.
--- @param base_cfg table|nil
--- @param override_cfg table|nil
--- @return table
function M.merge(base_cfg, override_cfg)
  local merged = {}

  if type(base_cfg) == "table" then
    for k, v in pairs(base_cfg) do
      merged[k] = v
    end
  end

  if type(override_cfg) ~= "table" then return merged end

  for k, v in pairs(override_cfg) do
    if k == "markers" then
      if type(v) == "table" then
        markers.append_markers(merged, v)
      end
    elseif k == "tile" and type(v) == "table" then
      local tile = {}
      if type(merged.tile) == "table" then
        for tk, tv in pairs(merged.tile) do tile[tk] = tv end
      end
      for tk, tv in pairs(v) do tile[tk] = tv end
      merged.tile = tile
    else
      merged[k] = v
    end
  end

  return merged
end

--- Build config from YAML metadata map.
--- @param mm table|nil
--- @return table|nil, string|nil
function M.from_meta(mm)
  if type(mm) ~= "table" then return {} end
  local cfg = {}

  if mm.center then cfg.center = meta.list_to_nums(mm.center) end
  if mm.zoom then cfg.zoom = meta.to_num(mm.zoom) end
  if mm.height then cfg.height = meta.to_str(mm.height) end
  if mm.width then cfg.width = meta.to_str(mm.width) end

  if mm.tile and type(mm.tile) == "table" then
    cfg.tile = {}
    for k, v in pairs(mm.tile) do
      local s = meta.to_str(v)
      cfg.tile[k] = tonumber(s) or s
    end
  end

  if mm.markers and type(mm.markers) == "table" then
    local out = {}
    for _, m_meta in ipairs(mm.markers) do
      local m = markers.from_meta_entry(m_meta)
      if m ~= nil then table.insert(out, m) end
    end
    markers.append_markers(cfg, out)
  end

  local markers_file = meta.to_nonempty_str(mm["markers-file"]) or meta.to_nonempty_str(mm.markers_file)
  local markers_sep = meta.to_nonempty_str(mm["markers-sep"]) or meta.to_nonempty_str(mm.markers_sep)
  local markers_err = markers.append_from_file(cfg, markers_file, markers_sep)
  if markers_err then return nil, markers_err end

  for k, v in pairs(mm) do
    if not M.SPECIAL[k] then
      local b = meta.to_bool(v)
      if b ~= nil then
        cfg[k] = b
      else
        cfg[k] = meta.to_num(v) or meta.to_str(v)
      end
    end
  end

  return cfg, nil
end

--- Build config from inline shortcode kwargs.
--- @param kwargs table
--- @return table|nil, string|nil
function M.from_kwargs(kwargs)
  local cfg = {}

  local function kstr(k)
    return meta.to_nonempty_str(kwargs[k])
  end

  local c_str = kstr("center")
  if c_str then cfg.center = markers.parse_coord_str(c_str) end

  local z = kstr("zoom")
  if z then cfg.zoom = tonumber(z) end

  local h = kstr("height")
  if h then cfg.height = h end

  local w = kstr("width")
  if w then cfg.width = w end

  local t_str = kstr("tile")
  if t_str then
    local ok, parsed = pcall(pandoc.json.decode, t_str)
    if ok and type(parsed) == "table" then
      local tile_obj = nil
      if parsed.url ~= nil then
        tile_obj = parsed
      elseif #parsed == 1 and type(parsed[1]) == "table" and parsed[1].url ~= nil then
        tile_obj = parsed[1]
      end

      if tile_obj ~= nil then
        cfg.tile = {}
        for k, v in pairs(tile_obj) do
          if type(v) == "string" then
            cfg.tile[k] = tonumber(v) or v
          else
            cfg.tile[k] = v
          end
        end
      end
    end
  end

  local m_str = kstr("markers")
  if m_str then
    local ok, parsed = pcall(pandoc.json.decode, m_str)
    if ok and type(parsed) == "table" then
      local out = {}
      for _, mj in ipairs(parsed) do
        local m = markers.from_json_entry(mj)
        if m ~= nil then table.insert(out, m) end
      end
      markers.append_markers(cfg, out)
    end
  end

  local markers_file = kstr("markers-file") or kstr("markers_file")
  local markers_sep = kstr("markers-sep") or kstr("markers_sep")
  local markers_err = markers.append_from_file(cfg, markers_file, markers_sep)
  if markers_err then return nil, markers_err end

  for k, v in pairs(kwargs) do
    if not M.SPECIAL[k] then
      local s = meta.to_str(v)
      local s_lc = s and s:lower() or nil
      if s_lc == "true" then
        cfg[k] = true
      elseif s_lc == "false" then
        cfg[k] = false
      else
        cfg[k] = tonumber(s) or s
      end
    end
  end

  return cfg, nil
end

return M
