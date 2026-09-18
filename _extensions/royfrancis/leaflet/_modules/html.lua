--- @module html
--- HTML rendering and escaping helpers.

local M = {}

--- Escape user-facing text for inline HTML rendering.
--- @param s any
--- @return string
function M.escape(s)
  return (tostring(s)
    :gsub("&", "&amp;")
    :gsub("<", "&lt;")
    :gsub(">", "&gt;")
    :gsub('"', "&quot;"))
end

--- Build a styled inline warning for extension errors.
--- @param msg string
--- @return pandoc.RawInline
function M.error_inline(msg)
  return pandoc.RawInline("html",
    '<span style="color:#c00;background:#fee;border:1px solid #c00;'
      .. 'padding:2px 6px;border-radius:3px;font-family:monospace;font-size:.9em;">'
      .. '&#x26A0; leaflet: ' .. M.escape(msg) .. '</span>')
end

return M
