--- @module json
--- Lightweight JSON encoder for map config serialization.

local M = {}

--- Encode a Lua value as JSON.
--- @param val any
--- @return string
function M.encode(val)
  if val == nil then return "null" end
  if type(val) == "boolean" then return tostring(val) end
  if type(val) == "number" then
    if val == math.floor(val) and math.abs(val) < 1e15 then
      return string.format("%d", val)
    end
    return tostring(val)
  end
  if type(val) == "string" then
    val = val
      :gsub('\\', '\\\\')
      :gsub('"', '\\"')
      :gsub('\n', '\\n')
      :gsub('\r', '\\r')
      :gsub('\t', '\\t')
    return '"' .. val .. '"'
  end
  if type(val) == "table" then
    local n = #val
    local is_arr = n > 0
    if is_arr then
      for i = 1, n do
        if val[i] == nil then
          is_arr = false
          break
        end
      end
    end

    if is_arr then
      local parts = {}
      for _, v in ipairs(val) do
        table.insert(parts, M.encode(v))
      end
      return "[" .. table.concat(parts, ",") .. "]"
    end

    local keys = {}
    for k in pairs(val) do
      table.insert(keys, k)
    end
    table.sort(keys, function(a, b) return tostring(a) < tostring(b) end)

    local parts = {}
    for _, k in ipairs(keys) do
      table.insert(parts, M.encode(tostring(k)) .. ":" .. M.encode(val[k]))
    end
    return "{" .. table.concat(parts, ",") .. "}"
  end

  return "null"
end

return M
