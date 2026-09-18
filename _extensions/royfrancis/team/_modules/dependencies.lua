--- @module dependencies
--- HTML and revealjs dependency registration helpers.

local M = {}

local html_deps_added = false
local revealjs_deps_added = false

--- Register HTML dependencies once per document.
function M.add_html_once()
  if html_deps_added then return end
  html_deps_added = true

  quarto.doc.add_html_dependency({
    name = "team",
    stylesheets = { "team.css" },
  })
end

--- Register revealjs dependencies once per document.
function M.add_revealjs_once()
  if revealjs_deps_added then return end
  revealjs_deps_added = true

  M.add_html_once()

  quarto.doc.add_html_dependency({
    name = "team-revealjs",
    stylesheets = { "team-revealjs.css" },
  })
end

return M