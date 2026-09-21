--- @module utils
--- Shared helpers for parsing, coercion, escaping, and block conversion.

local M = {}

--- Check if a value is empty after string conversion.
--- @param value any
--- @return boolean
function M.is_empty(value)
  return value == nil or pandoc.utils.stringify(value) == ""
end

--- Escape LaTeX special characters in text.
--- @param s string
--- @return string
function M.escape_latex(s)
  local out = s
  out = out:gsub("\\", "\0BACKSLASH\0")
  out = out:gsub("%%", "\\%%")
  out = out:gsub("{", "\\{")
  out = out:gsub("}", "\\}")
  out = out:gsub("#", "\\#")
  out = out:gsub("&", "\\&")
  out = out:gsub("%$", "\\$")
  out = out:gsub("_", "\\_")
  out = out:gsub("~", "\\textasciitilde{}")
  out = out:gsub("%^", "\\textasciicircum{}")
  out = out:gsub("\0BACKSLASH\0", "\\textbackslash{}")
  return out
end

--- Escape Typst special characters in text.
--- @param s string
--- @return string
function M.escape_typst(s)
  local out = s
  out = out:gsub("\\", "\\\\")
  out = out:gsub("%[", "\\[")
  out = out:gsub("%]", "\\]")
  out = out:gsub("_", "\\_")
  out = out:gsub("%*", "\\*")
  out = out:gsub("#", "\\#")
  out = out:gsub("%$", "\\$")
  out = out:gsub("@", "\\@")
  out = out:gsub("<", "\\<")
  out = out:gsub(">", "\\>")
  return out
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

--- Convert a value into Pandoc Blocks.
--- @param value any
--- @return pandoc.Blocks
function M.content_to_blocks(value)
  if value == nil then
    return pandoc.Blocks({})
  end

  local value_type = pandoc.utils.type(value)
  if value_type == "Inlines" then
    return pandoc.Blocks({ pandoc.Plain(value) })
  end
  if value_type == "Blocks" then
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

--- Wrap top-level plaintext and paragraphs in bold styling.
--- @param blocks pandoc.Blocks
--- @return pandoc.Blocks
function M.make_bold(blocks)
  local result = pandoc.Blocks({})
  for _, block in ipairs(blocks) do
    if block.t == "Plain" then
      result:insert(pandoc.Plain({ pandoc.Strong(block.content) }))
    elseif block.t == "Para" then
      result:insert(pandoc.Para({ pandoc.Strong(block.content) }))
    else
      result:insert(block)
    end
  end
  return result
end

--- Build a stable collapse id for one accordion item.
--- @param accordion_id string
--- @param item table
--- @param index integer
--- @return string
function M.build_collapse_id(accordion_id, item, index)
  if item.id ~= nil then
    local explicit = pandoc.utils.stringify(item.id)
    if explicit ~= "" then
      return "-" .. accordion_id .. "-" .. explicit
    end
  end

  local header = pandoc.utils.stringify(item.header or "")
  local body = pandoc.utils.stringify(item.body or "")
  local tail = (header:sub(-10) .. body:sub(-10)):lower():gsub("[^%w]", "")

  if tail == "" then
    tail = "item" .. tostring(index)
  else
    tail = tail .. "-" .. tostring(index)
  end

  return "-" .. accordion_id .. "-" .. tail
end

--- Build a missing-field error string for an item.
--- @param user_label string
--- @param index integer
--- @param header_content string
--- @param body_content string
--- @return string|nil
function M.item_validation_error(user_label, index, header_content, body_content)
  if header_content ~= "" and body_content ~= "" then
    return nil
  end

  local missing = {}
  if header_content == "" then table.insert(missing, "header") end
  if body_content == "" then table.insert(missing, "body") end
  local missing_str = table.concat(missing, "' and '")
  return string.format("'%s': Item %d is missing '%s'.", user_label, index, missing_str)
end

return M