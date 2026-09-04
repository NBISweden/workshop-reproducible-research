--- @module markers
--- Marker parsing, normalization, and file loading utilities.

local function load_module(path)
  return require(quarto.utils.resolve_path(path):gsub("%.lua$", ""))
end

local meta = load_module("_modules/metadata.lua")

local M = {}

--- Trim leading and trailing whitespace.
--- @param s any
--- @return string|nil
local function trim(s)
  if s == nil then return nil end
  return (tostring(s):gsub("^%s+", ""):gsub("%s+$", ""))
end

--- Parse a coordinate string in [lat, lon] form.
--- @param s string|nil
--- @return table|nil
function M.parse_coord_str(s)
  if s == nil then return nil end
  local inner = s:match("^%s*%[(.-)%]%s*$")
  if inner == nil then return nil end

  local nums = {}
  for part in inner:gmatch("[^,]+") do
    local n = tonumber(part:match("^%s*(.-)%s*$"))
    if n == nil then return nil end
    table.insert(nums, n)
  end

  return #nums == 2 and nums or nil
end

--- Append marker entries to cfg.markers while preserving existing values.
--- @param target table
--- @param markers table
function M.append_markers(target, markers)
  if type(markers) ~= "table" then return end
  if target.markers == nil then target.markers = {} end
  for _, marker in ipairs(markers) do
    table.insert(target.markers, marker)
  end
end

--- Extract marker coordinates from multiple supported marker schemas.
--- @param marker table
--- @return number|nil, number|nil
function M.marker_coords(marker)
  if type(marker) ~= "table" then return nil end

  local lat = tonumber(marker.lat)
  local lon = tonumber(marker.lon)
  if lat ~= nil and lon ~= nil then
    return lat, lon
  end

  if type(marker.position) == "table" and #marker.position == 2 then
    lat = tonumber(marker.position[1])
    lon = tonumber(marker.position[2])
    if lat ~= nil and lon ~= nil then
      return lat, lon
    end
  end

  return nil
end

--- Build a marker table from metadata map entry.
--- @param m_meta table
--- @return table|nil
function M.from_meta_entry(m_meta)
  if type(m_meta) ~= "table" then return nil end

  local m = {}
  local lat = meta.to_num(m_meta.lat or m_meta.latitude)
  local lon = meta.to_num(m_meta.lon or m_meta.lng or m_meta.long or m_meta.longitude)

  if lat ~= nil and lon ~= nil then
    m.lat = lat
    m.lon = lon
  elseif m_meta.position then
    local pos = meta.list_to_nums(m_meta.position)
    if pos and #pos == 2 then
      m.lat = pos[1]
      m.lon = pos[2]
    end
  end

  if m_meta.popup then m.popup = meta.to_html(m_meta.popup) end
  if m_meta.tooltip then m.tooltip = meta.to_inline_html(m_meta.tooltip) end

  if m_meta.icon then
    local iv = m_meta.icon
    if type(iv) == "table" and iv.url then
      m.icon = {
        url = meta.to_str(iv.url),
        size = iv.size and meta.list_to_nums(iv.size),
        anchor = iv.anchor and meta.list_to_nums(iv.anchor),
      }
    else
      m.icon = meta.to_str(iv)
    end
  end

  if m_meta["icon-color"] then m.icon_color = meta.to_str(m_meta["icon-color"]) end
  if m_meta["icon-size"] then m.icon_size = meta.to_num(m_meta["icon-size"]) end
  if m_meta["icon-anchor"] then m.icon_anchor = meta.list_to_nums(m_meta["icon-anchor"]) end

  return m
end

--- Build a marker table from decoded JSON entry.
--- @param mj table
--- @return table|nil
function M.from_json_entry(mj)
  if type(mj) ~= "table" then return nil end

  local m = {}
  local lat = tonumber(mj.lat or mj.latitude)
  local lon = tonumber(mj.lon or mj.lng or mj.long or mj.longitude)

  if lat ~= nil and lon ~= nil then
    m.lat = lat
    m.lon = lon
  elseif type(mj.position) == "table" and #mj.position == 2 then
    m.lat = tonumber(mj.position[1])
    m.lon = tonumber(mj.position[2])
  end

  m.popup = mj.popup
  m.tooltip = mj.tooltip
  m.icon = mj.icon
  m.icon_color = mj["icon-color"]
  m.icon_size = mj["icon-size"]
  m.icon_anchor = mj["icon-anchor"]

  return m
end

--- Compute center from marker bounding box midpoint.
--- @param markers table
--- @return table|nil
function M.center_from_markers(markers)
  if type(markers) ~= "table" then return nil end

  local count = 0
  local lat_min, lat_max = math.huge, -math.huge
  local lon_min, lon_max = math.huge, -math.huge

  for _, marker in ipairs(markers) do
    local lat, lon = M.marker_coords(marker)
    if lat ~= nil and lon ~= nil then
      count = count + 1
      if lat < lat_min then lat_min = lat end
      if lat > lat_max then lat_max = lat end
      if lon < lon_min then lon_min = lon end
      if lon > lon_max then lon_max = lon end
    end
  end

  if count == 0 then return nil end
  return { (lat_min + lat_max) / 2, (lon_min + lon_max) / 2 }
end

--- Split a delimited text row with optional CSV quoting support.
--- @param line string
--- @param separator string|nil
--- @return table
local function split_delimited_line(line, separator)
  local cells = {}
  line = (line or ""):gsub("\r$", "")
  local sep = separator or "\t"

  if sep == "\\t" then sep = "\t" end
  if sep == "" then
    table.insert(cells, trim(line) or "")
    return cells
  end

  if sep == "," then
    local cell = ""
    local in_quotes = false
    local i = 1

    while i <= #line do
      local ch = line:sub(i, i)
      if ch == '"' then
        if in_quotes and line:sub(i + 1, i + 1) == '"' then
          cell = cell .. '"'
          i = i + 1
        else
          in_quotes = not in_quotes
        end
      elseif ch == "," and not in_quotes then
        table.insert(cells, trim(cell) or "")
        cell = ""
      else
        cell = cell .. ch
      end
      i = i + 1
    end

    table.insert(cells, trim(cell) or "")
    return cells
  end

  local from = 1
  while true do
    local start_idx, end_idx = line:find(sep, from, true)
    if start_idx == nil then
      table.insert(cells, trim(line:sub(from)) or "")
      break
    end
    table.insert(cells, trim(line:sub(from, start_idx - 1)) or "")
    from = end_idx + 1
  end

  return cells
end

--- Infer separator from markers file extension.
--- @param path string|nil
--- @return string
local function default_separator_for_path(path)
  if path == nil then return "\t" end
  local p = tostring(path):lower()
  if p:match("%.csv$") then return "," end
  if p:match("%.tsv$") then return "\t" end
  return "\t"
end

--- Normalize delimiter header names for flexible matching.
--- @param key string
--- @return string
local function normalize_header_key(key)
  return ((trim(key) or ""):lower():gsub("[%s_%-]+", ""))
end

--- Resolve marker file path relative to source document.
--- @param path string
--- @return string
local function resolve_doc_path(path)
  if path == nil or path == "" then return path end
  if path:match("^/") or path:match("^%a:[/\\]") then return path end

  local input_files = PANDOC_STATE and PANDOC_STATE.input_files or nil
  local input_file = input_files and input_files[1] or nil
  if input_file == nil or input_file == "" then return path end

  local doc_dir = pandoc.path.directory(input_file)
  if doc_dir == nil or doc_dir == "" then return path end
  return pandoc.path.join({ doc_dir, path })
end

--- Build marker from normalized row values.
--- @param row table
--- @return table|nil
local function marker_from_row(row)
  local lat = tonumber(row.lat or row.latitude)
  local lon = tonumber(row.lon or row.lng or row.long or row.longitude)
  if lat == nil or lon == nil then return nil end

  local marker = { lat = lat, lon = lon }
  if row.popup and row.popup ~= "" then marker.popup = row.popup end
  if row.tooltip and row.tooltip ~= "" then marker.tooltip = row.tooltip end

  local icon_url = row.iconurl
  local icon = row.icon
  local icon_anchor = nil
  if row.iconanchor and row.iconanchor ~= "" then
    local x, y = tostring(row.iconanchor):match("(%d+),%s*(%d+)")
    if x and y then icon_anchor = { tonumber(x), tonumber(y) } end
  end
  local icon_size_pair = row.iconsize and M.parse_coord_str(row.iconsize) or nil

  if icon_url and icon_url ~= "" then
    marker.icon = { url = icon_url }
    if icon_size_pair then marker.icon.size = icon_size_pair end
    if icon_anchor then marker.icon.anchor = icon_anchor end
  elseif icon and icon ~= "" then
    marker.icon = icon
    if row.iconcolor and row.iconcolor ~= "" then marker.icon_color = row.iconcolor end
    if row.iconsize and row.iconsize ~= "" then marker.icon_size = tonumber(row.iconsize) end
    if icon_anchor then marker.icon_anchor = icon_anchor end
  end

  return marker
end

--- Parse a markers file into marker definitions.
--- @param markers_path string
--- @param separator string|nil
--- @return table|nil, string|nil
function M.from_file(markers_path, separator)
  local sep = separator or default_separator_for_path(markers_path)
  if sep == "\\t" then sep = "\t" end

  local resolved_path = resolve_doc_path(markers_path)
  local handle = nil
  local candidate_paths = { resolved_path }
  if resolved_path ~= markers_path then table.insert(candidate_paths, markers_path) end

  for _, candidate_path in ipairs(candidate_paths) do
    handle = io.open(candidate_path, "r")
    if handle ~= nil then break end
  end
  if handle == nil then
    return nil, string.format("could not read markers file '%s'", markers_path)
  end

  local header_keys = nil
  local markers = {}

  for line in handle:lines() do
    local trimmed = trim(line)
    if trimmed ~= nil and trimmed ~= "" and not trimmed:match("^#") then
      local cells = split_delimited_line(line, sep)
      if header_keys == nil then
        header_keys = {}
        for _, cell in ipairs(cells) do
          table.insert(header_keys, normalize_header_key(cell))
        end
      else
        local row = {}
        for idx, key in ipairs(header_keys) do
          row[key] = cells[idx] or ""
        end
        local marker = marker_from_row(row)
        if marker ~= nil then table.insert(markers, marker) end
      end
    end
  end

  handle:close()

  if header_keys == nil then
    return nil, string.format("markers file '%s' is empty", markers_path)
  end

  return markers
end

--- Load markers from file and append to config.
--- @param cfg table
--- @param markers_file string|nil
--- @param markers_sep string|nil
--- @return string|nil
function M.append_from_file(cfg, markers_file, markers_sep)
  if not markers_file then return nil end
  local markers, err = M.from_file(markers_file, markers_sep)
  if err then return err end
  M.append_markers(cfg, markers)
  return nil
end

return M
