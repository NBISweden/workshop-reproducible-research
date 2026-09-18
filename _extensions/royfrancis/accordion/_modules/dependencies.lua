--- @module dependencies
--- HTML and document dependency registration helpers.

local M = {}

local html_deps_added = false
local revealjs_deps_added = false
local latex_deps_added = false

--- Register HTML dependencies once per document.
function M.add_html_once()
  if html_deps_added then return end
  html_deps_added = true

  quarto.doc.add_html_dependency({
    name = "accordion",
    stylesheets = { "accordion.css" },
  })
end

--- Register revealjs dependencies once per document.
function M.add_revealjs_once()
  if revealjs_deps_added then return end
  revealjs_deps_added = true

  quarto.doc.add_html_dependency({
    name = "accordion",
    stylesheets = { "accordion.css" },
  })

  quarto.doc.add_html_dependency({
    name = "accordion-revealjs",
    stylesheets = { "accordion-revealjs.css" },
    scripts = { "accordion-revealjs.js" },
  })
end

--- Register LaTeX dependencies once per document.
function M.add_latex_once()
  if latex_deps_added then return end
  latex_deps_added = true

  quarto.doc.include_text(
    "in-header",
    "\\usepackage{tcolorbox}\n"
      .. "\\usepackage{xcolor}\n"
      .. "\\definecolor{accordionerror}{HTML}{721c24}"
  )
end

return M