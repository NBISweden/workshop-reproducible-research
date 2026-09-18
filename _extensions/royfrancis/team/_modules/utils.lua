--- @module utils
--- Shared helpers for parsing, coercion, escaping, and content conversion.

local M = {}

--- Check if a value is empty after string conversion.
--- @param value any
--- @return boolean
function M.is_empty(value)
  return value == nil or pandoc.utils.stringify(value) == ""
end

--- Strip one pair of surrounding single or double quotes.
--- @param s string|nil
--- @return string|nil
function M.strip_surrounding_quotes(s)
  if s == nil or #s < 2 then return s end

  local first = s:sub(1, 1)
  local last = s:sub(-1)
  if (first == '"' and last == '"') or (first == "'" and last == "'") then
    return s:sub(2, -2)
  end

  return s
end

--- Read a shortcode kwarg as an unquoted string.
--- @param kwargs table
--- @param key string
--- @return string
function M.get_kwarg(kwargs, key)
  if kwargs[key] == nil then return "" end
  return M.strip_surrounding_quotes(pandoc.utils.stringify(kwargs[key])) or ""
end

--- Validate shortcode labels.
--- @param label string
--- @return boolean
function M.is_valid_label(label)
  return label:find("^[%w%-_]+$") ~= nil
end

--- Produce a short deterministic hex string from content for auto-labelling.
--- Uses djb2 algorithm. Output is always 8 hex characters.
--- @param s string
--- @return string
function M.hash_string(s)
  local hash = 5381
  for i = 1, #s do
    hash = (hash * 33 + string.byte(s, i)) % 0x100000000
  end
  return string.format("%08x", hash)
end

--- Escape text for safe use in HTML attributes.
--- @param s string
--- @return string
function M.escape_html_attr(s)
  local out = s
  out = out:gsub("&", "&amp;")
  out = out:gsub('"', "&quot;")
  out = out:gsub("'", "&#39;")
  out = out:gsub("<", "&lt;")
  out = out:gsub(">", "&gt;")
  return out
end

--- Convert a value into Pandoc Blocks.
--- Strings are parsed as markdown so inline markdown and raw HTML are preserved.
--- @param value any
--- @return pandoc.Blocks
function M.content_to_blocks(value)
  if value == nil then
    return pandoc.Blocks({})
  end

  local value_type = pandoc.utils.type(value)
  if value_type == "Inlines" or value_type == "MetaInlines" then
    return pandoc.Blocks({ pandoc.Plain(value) })
  end
  if value_type == "Blocks" or value_type == "MetaBlocks" then
    return pandoc.Blocks(value)
  end
  if type(value) == "string" then
    local doc = pandoc.read(value, "markdown")
    if #doc.blocks == 1 and doc.blocks[1].t == "Para" then
      return pandoc.Blocks({ pandoc.Plain(doc.blocks[1].content) })
    end
    return doc.blocks
  end

  return pandoc.Blocks({
    pandoc.Plain({ pandoc.Str(pandoc.utils.stringify(value)) }),
  })
end

--- Convert a value into Pandoc Inlines.
--- @param value any
--- @return pandoc.Inlines
function M.content_to_inlines(value)
  if value == nil then
    return pandoc.Inlines({})
  end

  local value_type = pandoc.utils.type(value)
  if value_type == "Inlines" or value_type == "MetaInlines" then
    return pandoc.Inlines(value)
  end
  if value_type == "Blocks" or value_type == "MetaBlocks" then
    for _, block in ipairs(value) do
      if block.t == "Plain" or block.t == "Para" then
        return pandoc.Inlines(block.content)
      end
    end
    return pandoc.Inlines({ pandoc.Str(pandoc.utils.stringify(value)) })
  end
  if type(value) == "string" then
    local blocks = M.content_to_blocks(value)
    if #blocks > 0 and (blocks[1].t == "Plain" or blocks[1].t == "Para") then
      return pandoc.Inlines(blocks[1].content)
    end
    return pandoc.Inlines({ pandoc.Str(pandoc.utils.stringify(value)) })
  end

  return pandoc.Inlines({ pandoc.Str(pandoc.utils.stringify(value)) })
end

--- Build a stable DOM id for one team item.
--- If an explicit id is provided, that value is used as-is after string conversion.
--- @param team_id string
--- @param item table
--- @param index integer
--- @return string
function M.build_item_id(team_id, item, index)
  if item.id ~= nil then
    local explicit = pandoc.utils.stringify(item.id)
    if explicit ~= "" then
      return team_id .. "-" .. explicit
    end
  end

  local name = pandoc.utils.stringify(item.name or "")
  local description = pandoc.utils.stringify(item.description or "")
  local tail = (name .. description):lower():gsub("[^%w]", "")

  if tail == "" then
    tail = "member"
  else
    tail = tail:sub(1, 20)
  end

  return string.format("%s-%s-%d-%s", team_id, tail, index, M.hash_string(name .. description):sub(1, 6))
end

return M