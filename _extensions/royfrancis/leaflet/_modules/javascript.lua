--- @module javascript
--- JavaScript assembly for Leaflet map rendering.

local function load_module(path)
  return require(quarto.utils.resolve_path(path):gsub("%.lua$", ""))
end

local json = load_module("_modules/json.lua")
local markers = load_module("_modules/markers.lua")
local config = load_module("_modules/config.lua")

local M = {}

--- Build the JavaScript snippet for a map instance.
--- @param map_id string
--- @param map_var string
--- @param cfg table
--- @return string
function M.build(map_id, map_var, cfg)
  local out = {}

  table.insert(out,
    'if (typeof L.Icon.Default.imagePath !== "string") {\n'
      .. '  var _ql = document.querySelector(\'link[href*="quarto-leaflet"][href$="leaflet.css"]\');\n'
      .. '  if (_ql) { L.Icon.Default.imagePath = _ql.href.replace(/leaflet\\.css$/, "") + "images/"; }\n'
      .. '}\n')

  table.insert(out,
    'window.__quartoLeafletMaps = window.__quartoLeafletMaps || [];\n'
      .. 'window.__quartoLeafletWireReveal = window.__quartoLeafletWireReveal || function() {\n'
      .. '  if (window.__quartoLeafletRevealFix || typeof Reveal === "undefined") { return; }\n'
      .. '  window.__quartoLeafletRevealFix = true;\n'
      .. '  var _qlRefresh = function() {\n'
      .. '    setTimeout(function() {\n'
      .. '      window.__quartoLeafletMaps.forEach(function(m) {\n'
      .. '        if (m && typeof m.invalidateSize === "function") { m.invalidateSize(true); }\n'
      .. '      });\n'
      .. '    }, 60);\n'
      .. '  };\n'
      .. '  Reveal.on("ready", _qlRefresh);\n'
      .. '  Reveal.on("slidechanged", _qlRefresh);\n'
      .. '  Reveal.on("resize", _qlRefresh);\n'
      .. '  _qlRefresh();\n'
      .. '};\n'
      .. 'window.__quartoLeafletWireReveal();\n'
      .. 'if (document.readyState === "loading") {\n'
      .. '  document.addEventListener("DOMContentLoaded", window.__quartoLeafletWireReveal, {once: true});\n'
      .. '} else {\n'
      .. '  setTimeout(window.__quartoLeafletWireReveal, 0);\n'
      .. '}\n')

  local map_opts = { center = cfg.center, zoom = cfg.zoom }
  for k, v in pairs(cfg) do
    if not config.SPECIAL[k] then map_opts[k] = v end
  end

  table.insert(out, string.format("var %s = L.map('%s', %s);\n", map_var, map_id, json.encode(map_opts)))

  local tile_url = "https://tile.openstreetmap.org/{z}/{x}/{y}.png"
  local tile_opts = {
    attribution = '&copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> contributors'
  }

  if cfg.tile then
    if cfg.tile.url then tile_url = cfg.tile.url end
    if cfg.tile.attribution then tile_opts.attribution = cfg.tile.attribution end
    for k, v in pairs(cfg.tile) do
      if k ~= "url" and k ~= "attribution" then tile_opts[k] = v end
    end
  end

  table.insert(out, string.format("L.tileLayer(%s, %s).addTo(%s);\n", json.encode(tile_url), json.encode(tile_opts), map_var))

  if cfg.markers then
    for _, m in ipairs(cfg.markers) do
      local lat, lon = markers.marker_coords(m)
      if lat ~= nil and lon ~= nil then
        local coords = { lat, lon }
        local popup_offset = nil
        local tooltip_offset = nil

        if m.icon then
          local ic = m.icon
          if type(ic) == "table" and ic.url then
            local opts = {
              iconUrl = ic.url,
              iconSize = ic.size or { 25, 41 },
              iconAnchor = m.icon_anchor or ic.anchor or { 12, 41 },
            }
            table.insert(out, string.format("var _icon = L.icon(%s);\n", json.encode(opts)))
            table.insert(out, string.format("var _m = L.marker(%s, {icon: _icon}).addTo(%s);\n", json.encode(coords), map_var))
          else
            local sz = m.icon_size or 14
            local color = m.icon_color or "currentColor"
            local pin_padding = 4
            local pin_head_ratio = 0.6
            local pin_size = math.max(math.ceil((sz + 2 * pin_padding) / pin_head_ratio), sz + 10)
            local ax = m.icon_anchor and m.icon_anchor[1] or math.floor(pin_size / 2)
            local ay = m.icon_anchor and m.icon_anchor[2] or pin_size
            local html = string.format(
              '<span class="quarto-leaflet-icon-pin" style="--ql-icon-size:%dpx;--ql-pin-size:%dpx;--ql-pin-color:%s;"><i class="fa-solid fa-location-pin quarto-leaflet-icon-pin-bg" aria-hidden="true"></i><i class="%s quarto-leaflet-icon-glyph" aria-hidden="true"></i></span>',
              sz,
              pin_size,
              color,
              ic
            )
            table.insert(out, string.format(
              'var _icon = L.divIcon({html: %s, className: "leaflet-div-icon-custom", iconSize: [%d, %d], iconAnchor: [%d, %d]});\n',
              json.encode(html), pin_size, pin_size, ax, ay
            ))
            table.insert(out, string.format("var _m = L.marker(%s, {icon: _icon}).addTo(%s);\n", json.encode(coords), map_var))
            popup_offset = { 0, -pin_size }
            tooltip_offset = { math.ceil(pin_size / 2), -math.ceil(pin_size / 2) }
          end
        else
          table.insert(out, string.format("var _m = L.marker(%s).addTo(%s);\n", json.encode(coords), map_var))
        end

        if m.popup then
          if popup_offset then
            table.insert(out, string.format("_m.bindPopup(%s, {offset: %s});\n", json.encode(m.popup), json.encode(popup_offset)))
          else
            table.insert(out, string.format("_m.bindPopup(%s);\n", json.encode(m.popup)))
          end
        end

        if m.tooltip then
          if tooltip_offset then
            table.insert(out, string.format("_m.bindTooltip(%s, {offset: %s});\n", json.encode(m.tooltip), json.encode(tooltip_offset)))
          else
            table.insert(out, string.format("_m.bindTooltip(%s);\n", json.encode(m.tooltip)))
          end
        end
      end
    end
  end

  table.insert(out, string.format("window.__quartoLeafletMaps.push(%s);\n", map_var))

  return table.concat(out)
end

return M
