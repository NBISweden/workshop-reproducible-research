--- @module leaflet
--- @license MIT
--- @author Roy Francis
--- Quarto shortcode extension entrypoint for Leaflet maps.

local function load_module(path)
  return require(quarto.utils.resolve_path(path):gsub("%.lua$", ""))
end

local cfg_mod = load_module("_modules/config.lua")
local js_mod = load_module("_modules/javascript.lua")
local html_mod = load_module("_modules/html.lua")
local markers_mod = load_module("_modules/markers.lua")
local resources_mod = load_module("_modules/resources.lua")

local map_counter = 0

--- Parse shortcode arguments into extension config.
--- @param args table
--- @param kwargs table
--- @param meta table
--- @return table|nil, pandoc.RawInline|nil
local function parse_cfg(args, kwargs, meta)
  if #args == 1 and next(kwargs) == nil then
    local label = pandoc.utils.stringify(args[1])
    local leaflet_meta = meta.leaflet
    if leaflet_meta == nil then
      return nil, html_mod.error_inline("no 'leaflet' key in document metadata")
    end

    local map_meta = leaflet_meta[label]
    if map_meta == nil then
      return nil, html_mod.error_inline(string.format("map '%s' not found in metadata", label))
    end

    local defaults_cfg, defaults_err = cfg_mod.from_meta(cfg_mod.build_leaflet_defaults_meta(leaflet_meta))
    if defaults_err then return nil, html_mod.error_inline(defaults_err) end

    local map_cfg, cfg_err = cfg_mod.from_meta(map_meta)
    if cfg_err then return nil, html_mod.error_inline(cfg_err) end

    return cfg_mod.merge(defaults_cfg, map_cfg), nil
  end

  local cfg, cfg_err = cfg_mod.from_kwargs(kwargs)
  if cfg_err then
    return nil, html_mod.error_inline(cfg_err)
  end
  return cfg, nil
end

--- Render a Leaflet shortcode instance.
--- @param args table
--- @param kwargs table
--- @param meta table
--- @return pandoc.RawInline|pandoc.RawBlock
local function render_leaflet(args, kwargs, meta)
  local cfg, cfg_error = parse_cfg(args, kwargs, meta)
  if cfg_error then return cfg_error end

  if cfg.center == nil then
    cfg.center = markers_mod.center_from_markers(cfg.markers)
  end
  if cfg.zoom == nil then
    cfg.zoom = 13
  end
  if cfg.center == nil then
    return html_mod.error_inline("'center' is required")
  end

  if not (quarto.doc.is_format("html") or quarto.doc.is_format("revealjs")) then
    return pandoc.RawBlock("latex", "")
  end

  resources_mod.add_once()

  map_counter = map_counter + 1
  local map_id = "quarto-leaflet-map-" .. map_counter
  local map_var = "quartoLeafletMap" .. map_counter

  local height = cfg.height or "400px"
  local width = cfg.width or "100%"
  local div = string.format(
    '<div id="%s" class="quarto-leaflet-map" style="height: %s; width: %s;"></div>',
    map_id,
    height,
    width
  )

  local ok_js, js_result = pcall(js_mod.build, map_id, map_var, cfg)
  if not ok_js then
    return html_mod.error_inline("render error: " .. tostring(js_result))
  end

  return pandoc.RawInline("html", div .. "\n" .. string.format("<script>\n%s</script>", js_result))
end

return {
  leaflet = render_leaflet,
}
