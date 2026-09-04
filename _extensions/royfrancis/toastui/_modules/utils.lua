--- Utility helpers for metadata conversion, shortcode args, and JSON serialization.
--- @module toastui._modules.utils

local M = {}

--- Convert a Pandoc value to plain text.
--- @param val any
--- @return string|nil
function M.stringify(val)
  if val == nil then return nil end
  return pandoc.utils.stringify(val)
end

--- Read a shortcode kwarg and normalize empty values to nil.
--- @param kwargs table
--- @param key string
--- @return string|nil
function M.get_kwarg(kwargs, key)
  local val = kwargs[key]
  if val == nil then return nil end
  local s = M.stringify(val)
  if s == nil or s == "" then return nil end
  s = s:gsub('^"(.*)"$', '%1'):gsub("^'(.*)'$", '%1')
  if s == "" then return nil end
  return s
end

--- Escape text for HTML attribute context.
--- @param s any
--- @return string
function M.escape_html_attr(s)
  if s == nil then return "" end
  s = tostring(s)
  s = s:gsub("&", "&amp;")
  s = s:gsub('"', "&quot;")
  s = s:gsub("'", "&#39;")
  s = s:gsub("<", "&lt;")
  s = s:gsub(">", "&gt;")
  return s
end

--- Convert Pandoc metadata values to plain Lua values.
--- @param val any
--- @return any
function M.meta_to_lua(val)
  if val == nil then return nil end
  local t = pandoc.utils.type(val)
  if t == "Inlines" or t == "Blocks" or t == "MetaInlines" or t == "MetaBlocks" then
    local s = pandoc.utils.stringify(val)
    if s == "true" then return true end
    if s == "false" then return false end
    local num = tonumber(s)
    if num then return num end
    return s
  elseif t == "string" or t == "boolean" or t == "number" then
    return val
  elseif t == "MetaMap" or t == "table" then
    local result = {}
    for k, v in pairs(val) do
      result[k] = M.meta_to_lua(v)
    end
    return result
  elseif t == "MetaList" or t == "List" then
    local result = {}
    for i, v in ipairs(val) do
      result[i] = M.meta_to_lua(v)
    end
    return result
  else
    return pandoc.utils.stringify(val)
  end
end

--- Serialize Lua value into JSON-safe JavaScript literal text.
--- @param val any
--- @return string
function M.to_json(val)
  if type(val) == "table" then
    return quarto.json.encode(val)
  elseif type(val) == "string" then
    if val == "true" then return "true" end
    if val == "false" then return "false" end
    local num = tonumber(val)
    if num then return tostring(num) end
    return quarto.json.encode(val)
  elseif type(val) == "boolean" then
    return val and "true" or "false"
  elseif type(val) == "number" then
    return tostring(val)
  elseif val == nil then
    return "null"
  else
    return quarto.json.encode(tostring(val))
  end
end

--- Resolve a path relative to input document directory when needed.
--- @param filepath string
--- @param doc_dir string|nil
--- @return string
function M.resolve_input_relative_path(filepath, doc_dir)
  if doc_dir and not filepath:match("^/") then
    return doc_dir .. "/" .. filepath
  end
  return filepath
end

--- Normalize escaped separator sequences from YAML/shortcode.
--- @param sep string|nil
--- @return string
function M.normalize_separator(sep)
  local out = sep or "\t"
  return out:gsub("\\t", "\t"):gsub("\\n", "\n"):gsub("\\\\", "\\")
end

--- Convert numeric-like height into CSS px value.
--- @param height any
--- @return string
function M.normalize_height(height)
  local out = height or "600px"
  if tonumber(out) then
    out = out .. "px"
  end
  return out
end

--- Return true when current output supports interactive calendar rendering.
--- @return boolean
function M.is_supported_format()
  return quarto.doc.is_format("html:js") or quarto.doc.is_format("revealjs")
end

--- Parse a shortcode kwarg string value into a typed Lua value.
--- Tries JSON decode first (for objects/arrays), then boolean, then number, else string.
--- @param s string
--- @return any
function M.parse_kwarg_value(s)
  if s == nil then return nil end
  local first = s:sub(1, 1)
  if first == "{" or first == "[" then
    local ok, decoded = pcall(quarto.json.decode, s)
    if ok and decoded ~= nil then
      return decoded
    end
  end
  if s == "true" then return true end
  if s == "false" then return false end
  local num = tonumber(s)
  if num then return num end
  return s
end

return M
