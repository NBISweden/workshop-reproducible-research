--- @module resources
--- Resource registration for HTML/revealjs output.

local M = {}

local resources_added = false
local fontawesome_added = false

--- Ensure extension resources are added at most once per document.
function M.add_once()
  if resources_added then return end
  resources_added = true

  if not fontawesome_added then
    quarto.doc.include_text(
      "in-header",
      '<link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.5.2/css/all.min.css">'
    )
    fontawesome_added = true
  end

  quarto.doc.add_html_dependency({
    name = "quarto-leaflet",
    version = "1.0.0",
    resources = {
      { path = "assets/images/marker-icon.png", name = "images/marker-icon.png" },
      { path = "assets/images/marker-icon-2x.png", name = "images/marker-icon-2x.png" },
      { path = "assets/images/marker-shadow.png", name = "images/marker-shadow.png" },
      { path = "assets/images/layers.png", name = "images/layers.png" },
      { path = "assets/images/layers-2x.png", name = "images/layers-2x.png" },
    },
    stylesheets = { "assets/leaflet.css", "quarto-leaflet.css" },
    scripts = { { path = "assets/leaflet.js" } },
  })
end

return M
