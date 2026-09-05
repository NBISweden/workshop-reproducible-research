--- Dependency registration for TOAST UI assets.
--- @module toastui._modules.dependencies

local M = {}
local deps_added = false

--- Register JS/CSS dependencies exactly once per document render.
--- @return nil
function M.add_deps_once()
  if deps_added then return end
  deps_added = true
  quarto.doc.add_html_dependency({
    name = "toastui-calendar",
    version = "2.1.3",
    scripts = { "assets/toastui-calendar.min.js" },
    stylesheets = {
      "assets/toastui-calendar.min.css",
      "toastui.css",
    },
  })
end

return M
