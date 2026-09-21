--- @module toastui
--- Quarto shortcode extension entrypoint for TOAST UI Calendar.

local deps = require('./_modules/dependencies')
local utils = require('./_modules/utils')
local config = require('./_modules/config')
local events = require('./_modules/events')
local render = require('./_modules/render')

--- Render the toastui shortcode.
--- @param args table
--- @param kwargs table
--- @param meta table
--- @param raw_args any
--- @param context string
--- @return PandocRawBlock|PandocNull
local function render_shortcode(args, kwargs, meta, raw_args, context)
  if not utils.is_supported_format() then
    return pandoc.Null()
  end

  local ok, result = pcall(function()
    deps.add_deps_once()

    local doc_dir = config.get_document_dir()
    local label = config.get_label(args)
    local cfg = config.build_cfg(label, meta, kwargs)

    local opts = config.build_options(cfg)
    local calendars = config.build_calendars(cfg)
    local event_list, event_errors = events.build_events(cfg, doc_dir)
    event_errors = event_errors or {}
    local initial_date = config.initial_date(cfg)
    local time_format, time_format_error = config.time_format(cfg)
    if time_format_error then
      table.insert(event_errors, time_format_error)
    end
    local event_detail_items, event_detail_items_error = config.event_detail_items(cfg)
    if event_detail_items_error then
      table.insert(event_errors, event_detail_items_error)
    end
    local popup_detail_items, popup_detail_items_error = config.popup_detail_items(cfg)
    if popup_detail_items_error then
      table.insert(event_errors, popup_detail_items_error)
    end

    local show_nav = config.show_nav(cfg)
    local height = utils.normalize_height(cfg.height)
    local timegrid_height = cfg.timegridHeight or "200%"

    return {
      opts = opts,
      calendars = calendars,
      events = event_list,
      initial_date = initial_date,
      time_format = time_format,
      event_detail_items = event_detail_items,
      popup_detail_items = popup_detail_items,
      show_nav = show_nav,
      height = height,
      timegrid_height = timegrid_height,
      errors = event_errors or {},
    }
  end)

  if not ok then
    return render.render_error_block({ "toastui: " .. tostring(result) })
  end

  if type(result.errors) == "table" and #result.errors > 0 then
    return render.render_error_block(result.errors)
  end

  local render_ok, rendered = pcall(function()
    return render.render_calendar_block(
      result.opts,
      result.calendars,
      result.events,
      result.initial_date,
      result.show_nav,
      result.height,
      result.timegrid_height,
      result.time_format,
      result.event_detail_items,
      result.popup_detail_items
    )
  end)

  if not render_ok then
    return render.render_error_block({ "toastui: " .. tostring(rendered) })
  end

  return rendered
end

return {
  ["toastui"] = render_shortcode,
}
