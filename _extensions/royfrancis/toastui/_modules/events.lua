--- Event loading, parsing, and validation utilities.
--- @module toastui._modules.events

local utils = require('./utils')

local M = {}
local REQUIRED_EVENT_COLUMNS = { "title", "start", "end" }

--- Split a string using a plain-text separator.
--- @param str string
--- @param sep string
--- @return table
local function split_line(str, sep)
  local parts = {}
  local start = 1
  while true do
    local pos = str:find(sep, start, true)
    if pos then
      table.insert(parts, str:sub(start, pos - 1))
      start = pos + #sep
    else
      table.insert(parts, str:sub(start))
      break
    end
  end
  return parts
end

--- Validate required event fields for a normalized event list.
--- @param events table
--- @param source string
--- @return table
local function validate_event_list(events, source)
  local errors = {}

  if type(events) ~= "table" then
    table.insert(errors, "toastui: events from " .. source .. " are not a table")
    return errors
  end

  for idx, ev in ipairs(events) do
    if type(ev) ~= "table" then
      table.insert(errors, "toastui: event #" .. idx .. " from " .. source .. " is not an object")
    else
      if ev.title == nil or tostring(ev.title) == "" then
        table.insert(errors, "toastui: event #" .. idx .. " from " .. source .. " is missing required field 'title'")
      end
      if ev.start == nil or tostring(ev.start) == "" then
        table.insert(errors, "toastui: event #" .. idx .. " from " .. source .. " is missing required field 'start'")
      end
      if ev["end"] == nil or tostring(ev["end"]) == "" then
        table.insert(errors, "toastui: event #" .. idx .. " from " .. source .. " is missing required field 'end'")
      end
    end
  end

  return errors
end

--- Validate required header columns for delimited event files.
--- @param headers table
--- @param source string
--- @return table
local function validate_required_columns(headers, source)
  local present = {}
  local missing = {}

  for _, key in ipairs(headers) do
    if key and key ~= "" then
      present[key] = true
    end
  end

  for _, required in ipairs(REQUIRED_EVENT_COLUMNS) do
    if not present[required] then
      table.insert(missing, required)
    end
  end

  if #missing == 0 then
    return {}
  end

  return {
    "toastui: events file " .. source .. " is missing required columns: " .. table.concat(missing, ", "),
  }
end

--- Parse a delimited event file with a header row into EventObject-like entries.
--- @param filepath string
--- @param sep string
--- @return table|nil,table
local function parse_events_file(filepath, sep)
  local errors = {}
  local file = io.open(filepath, "r")
  if not file then
    table.insert(errors, "toastui: cannot open events file: " .. filepath)
    return nil, errors
  end

  local content = file:read("*a")
  file:close()

  if not content or content == "" then
    table.insert(errors, "toastui: events file is empty: " .. filepath)
    return nil, errors
  end

  local lines = {}
  for line in content:gmatch("[^\r\n]+") do
    table.insert(lines, line)
  end

  if #lines < 2 then
    table.insert(errors, "toastui: events file has no data rows: " .. filepath)
    return nil, errors
  end

  local headers = {}
  local header_parts = split_line(lines[1], sep)
  for _, col in ipairs(header_parts) do
    table.insert(headers, col:match("^%s*(.-)%s*$"))
  end

  local header_errors = validate_required_columns(headers, filepath)
  for _, err in ipairs(header_errors) do
    table.insert(errors, err)
  end
  if #errors > 0 then
    return nil, errors
  end

  local numeric_keys = {
    goingDuration = true, comingDuration = true,
  }

  local events = {}
  for i = 2, #lines do
    local event = {}
    local parts = split_line(lines[i], sep)
    for col_idx, val in ipairs(parts) do
      local trimmed = val:match("^%s*(.-)%s*$")
      if col_idx <= #headers then
        local key = headers[col_idx]
        if trimmed == "true" then
          event[key] = true
        elseif trimmed == "false" then
          event[key] = false
        elseif key == "attendees" then
          if trimmed ~= "" then
            event[key] = { trimmed }
          else
            event[key] = {}
          end
        elseif numeric_keys[key] then
          local num = tonumber(trimmed)
          if num then
            event[key] = num
          end
        elseif trimmed ~= "" then
          event[key] = trimmed
        end
      end
    end
    table.insert(events, event)
  end

  local event_errors = validate_event_list(events, filepath)
  for _, err in ipairs(event_errors) do
    table.insert(errors, err)
  end

  return events, errors
end

--- Build the effective events list from YAML inline data and/or file data.
--- @param cfg table
--- @param doc_dir string|nil
--- @return table|nil,table
function M.build_events(cfg, doc_dir)
  local events = nil
  local errors = {}

  if cfg.events then
    events = cfg.events
    for _, ev in ipairs(events) do
      if type(ev) == "table" and ev.attendees ~= nil and type(ev.attendees) ~= "table" then
        ev.attendees = { tostring(ev.attendees) }
      end
    end
    local metadata_errors = validate_event_list(events, "metadata")
    for _, err in ipairs(metadata_errors) do
      table.insert(errors, err)
    end
  end

  local filepath = cfg.file
  if filepath and filepath ~= "" then
    local sep = utils.normalize_separator(cfg["file-sep"])
    local full_path = utils.resolve_input_relative_path(filepath, doc_dir)
    local file_events, file_errors = parse_events_file(full_path, sep)
    events = file_events
    for _, err in ipairs(file_errors) do
      table.insert(errors, err)
    end
  end

  return events, errors
end

M.validate_event_list = validate_event_list
M.validate_required_columns = validate_required_columns
M.parse_events_file = parse_events_file

return M
