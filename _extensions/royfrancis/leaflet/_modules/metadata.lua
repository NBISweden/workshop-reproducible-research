--- @module metadata
--- Metadata coercion helpers for Pandoc Meta values.

local M = {}

--- Convert any metadata value to a plain string.
--- @param val any Metadata value
--- @return string|nil
function M.to_str(val)
  if val == nil then return nil end
  if type(val) == "string" then return val end
  if type(val) == "number" then return tostring(val) end
  if type(val) == "boolean" then return tostring(val) end
  return pandoc.utils.stringify(val)
end

--- Convert metadata (including inlines/blocks) to HTML.
--- @param val any Metadata value
--- @return string|nil
function M.to_html(val)
  if val == nil then return nil end
  if type(val) == "string" then return val end
  local ok, html = pcall(function()
    return pandoc.write(pandoc.Pandoc({ pandoc.Para(val) }), "html")
  end)
  if ok and html then return html end
  return pandoc.utils.stringify(val)
end

--- Convert metadata to inline HTML (no paragraph wrappers).
--- @param val any Metadata value
--- @return string|nil
function M.to_inline_html(val)
  if val == nil then return nil end
  if type(val) == "string" then return val end
  local ok, html = pcall(function()
    return pandoc.write(pandoc.Pandoc({ pandoc.Plain(val) }), "html")
  end)
  if ok and html then
    return tostring(html):gsub("%s+$", "")
  end
  return pandoc.utils.stringify(val)
end

--- Convert metadata to number.
--- @param val any Metadata value
--- @return number|nil
function M.to_num(val)
  local s = M.to_str(val)
  return s and tonumber(s)
end

--- Convert metadata to non-empty trimmed string.
--- @param val any Metadata value
--- @return string|nil
function M.to_nonempty_str(val)
  local s = M.to_str(val)
  if s ~= nil then
    s = tostring(s):gsub("^%s+", ""):gsub("%s+$", "")
  end
  if s == nil or s == "" then return nil end
  return s
end

--- Convert metadata to boolean.
--- @param val any Metadata value
--- @return boolean|nil
function M.to_bool(val)
  if type(val) == "boolean" then return val end
  local s = M.to_str(val)
  if s == "true" then return true end
  if s == "false" then return false end
  return nil
end

--- Convert a metadata list of scalar values into numbers.
--- @param val any Metadata list
--- @return table|nil
function M.list_to_nums(val)
  if type(val) ~= "table" then return nil end
  local result = {}
  for _, item in ipairs(val) do
    local n = M.to_num(item)
    if n == nil then return nil end
    table.insert(result, n)
  end
  return #result > 0 and result or nil
end

return M
