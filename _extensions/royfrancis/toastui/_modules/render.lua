--- Rendering utilities for calendar HTML and JavaScript payload.
--- @module toastui._modules.render

local utils = require('./utils')

local M = {}
local calendar_counter = 0
local ERROR_BOX_STYLE = "padding: 0.75rem 1rem; border: 1px solid #f5c6cb; border-radius: 6px; color: #721c24; background: #f8d7da;"

--- Generate a unique DOM id for each calendar instance.
--- @return string
local function next_calendar_id()
  calendar_counter = calendar_counter + 1
  return "toastui-calendar-" .. calendar_counter
end

--- Build the navigation toolbar HTML.
--- @param container_id string
--- @return table
local function nav_html(container_id)
  return {
    '<div class="toastui-calendar-wrapper">',
    '<div class="toastui-nav" id="' .. utils.escape_html_attr(container_id) .. '-nav">',
    '  <button class="toastui-nav-btn" data-action="prev">&#9664;</button>',
    '  <button class="toastui-nav-btn" data-action="today">Today</button>',
    '  <button class="toastui-nav-btn" data-action="next">&#9654;</button>',
    '  <span class="toastui-nav-title" id="' .. utils.escape_html_attr(container_id) .. '-title"></span>',
    '  <span class="toastui-nav-right">',
    '    <button class="toastui-nav-btn toastui-view-btn" data-view="month">Month</button>',
    '    <button class="toastui-nav-btn toastui-view-btn" data-view="week">Week</button>',
    '    <button class="toastui-nav-btn toastui-view-btn" data-view="day">Day</button>',
    '  </span>',
    '</div>',
  }
end

--- Build JavaScript code that updates title and binds nav actions.
--- @param container_id string
--- @return string
local function nav_js(container_id)
  return [[
  function updateTitle() {
    var date = cal.getDate().toDate();
    var viewName = cal.getViewName();
    var title = '';
    var months = ['January','February','March','April','May','June','July','August','September','October','November','December'];
    if (viewName === 'month') {
      title = months[date.getMonth()] + ' ' + date.getFullYear();
    } else if (viewName === 'week') {
      var start = cal.getDateRangeStart().toDate();
      var end_ = cal.getDateRangeEnd().toDate();
      title = months[start.getMonth()] + ' ' + start.getDate() + ' - ' + months[end_.getMonth()] + ' ' + end_.getDate() + ', ' + end_.getFullYear();
    } else {
      title = months[date.getMonth()] + ' ' + date.getDate() + ', ' + date.getFullYear();
    }
    var el = document.getElementById(]] .. utils.to_json(container_id .. "-title") .. [[);
    if (el) el.textContent = title;
    var navEl = document.getElementById(]] .. utils.to_json(container_id .. "-nav") .. [[);
    if (navEl) {
      var btns = navEl.querySelectorAll('.toastui-view-btn');
      for (var i = 0; i < btns.length; i++) {
        btns[i].classList.toggle('active', btns[i].getAttribute('data-view') === viewName);
      }
    }
  }
  var navEl = document.getElementById(]] .. utils.to_json(container_id .. "-nav") .. [[);
  if (navEl) {
    navEl.addEventListener('click', function(e) {
      var btn = e.target.closest('[data-action]');
      if (btn) {
        var action = btn.getAttribute('data-action');
        if (action === 'prev') cal.prev();
        else if (action === 'next') cal.next();
        else if (action === 'today') cal.today();
        updateTitle();
        return;
      }
      var viewBtn = e.target.closest('[data-view]');
      if (viewBtn) {
        cal.changeView(viewBtn.getAttribute('data-view'));
        updateTitle();
      }
    });
  }
  updateTitle();
]]
end

--- Render one or more error messages in the output document.
--- @param errors string|table
--- @return PandocRawBlock
function M.render_error_block(errors)
  local items = {}

  if type(errors) == "string" then
    table.insert(items, errors)
  elseif type(errors) == "table" then
    for _, message in ipairs(errors) do
      if message and tostring(message) ~= "" then
        table.insert(items, tostring(message))
      end
    end
  end

  if #items == 0 then
    table.insert(items, "toastui: an unexpected error occurred")
  end

  local function format_msg(raw)
    local escaped = utils.escape_html_attr(raw):gsub("&#10;", "<br>")
    return escaped:gsub("^toastui:", "<strong>toastui</strong>:")
  end

  if #items == 1 then
    return pandoc.RawBlock("html", '<div style="' .. ERROR_BOX_STYLE .. '">' .. format_msg(items[1]) .. '</div>')
  end

  local html = {
    '<div style="' .. ERROR_BOX_STYLE .. '">',
    '  <ul style="margin: 0; padding-left: 1.2rem;">',
  }

  for _, message in ipairs(items) do
    table.insert(html, '    <li>' .. format_msg(message) .. '</li>')
  end

  table.insert(html, '  </ul>')
  table.insert(html, '</div>')

  return pandoc.RawBlock("html", table.concat(html, "\n"))
end

--- Render a complete calendar widget as a raw HTML block.
--- @param opts table
--- @param calendars table|nil
--- @param events table|nil
--- @param initial_date string|nil
--- @param show_nav boolean
--- @param height string
--- @param timegrid_height string|nil
--- @return PandocRawBlock
function M.render_calendar_block(opts, calendars, events, initial_date, show_nav, height, timegrid_height)
  local container_id = next_calendar_id()
  local html_parts = {}

  local tg = timegrid_height or "200%"
  table.insert(html_parts, '<style>#' .. utils.escape_html_attr(container_id) .. ' .toastui-calendar-timegrid { height: ' .. utils.escape_html_attr(tg) .. '; min-height: unset; }</style>')

  if show_nav then
    local nav = nav_html(container_id)
    for _, line in ipairs(nav) do
      table.insert(html_parts, line)
    end
  else
    table.insert(html_parts, '<div class="toastui-calendar-wrapper">')
  end

  local detail_popup_class = ""
  if opts.useDetailPopup then
    detail_popup_class = " toastui-calendar-detail-popup-enabled"
  end
  table.insert(html_parts, '<div id="' .. utils.escape_html_attr(container_id) .. '" class="' .. detail_popup_class .. '" style="height: ' .. utils.escape_html_attr(height) .. ';"></div>')
  table.insert(html_parts, '</div>')

  table.insert(html_parts, '<script>')
  table.insert(html_parts, '(function() {')
  table.insert(html_parts, '  var opts = ' .. utils.to_json(opts) .. ';')

  if calendars then
    table.insert(html_parts, '  opts.calendars = ' .. utils.to_json(calendars) .. ';')
  end

  table.insert(html_parts, '  var cal = new tui.Calendar(document.getElementById(' .. utils.to_json(container_id) .. '), opts);')

  if events then
    -- Patch events so that detail-popup sections only appear for
    -- user-provided properties.  The library internally defaults
    -- state → "Busy" and attendees → [], both truthy, causing those
    -- sections to always render.  Setting explicit falsy values
    -- before createEvents() prevents those defaults from kicking in.
    table.insert(html_parts, '  var __ev = ' .. utils.to_json(events) .. ';')
    table.insert(html_parts, [[  __ev.forEach(function(e) {
    if (!e.state) e.state = '';
    if (!e.attendees || (Array.isArray(e.attendees) && e.attendees.length === 0)) e.attendees = null;
  });]])
    table.insert(html_parts, '  cal.createEvents(__ev);')
  end

  if initial_date and initial_date ~= "" then
    table.insert(html_parts, '  cal.setDate(new Date(' .. utils.to_json(initial_date) .. '));')
  end

  table.insert(html_parts, '  window.__quartoToastuiCalendars = window.__quartoToastuiCalendars || {};')
  table.insert(html_parts, '  window.__quartoToastuiCalendars[' .. utils.to_json(container_id) .. '] = cal;')
  table.insert(html_parts, '  window.__quartoToastuiLastCalendarId = ' .. utils.to_json(container_id) .. ';')

  -- Dismiss detail/form popup on click outside.
  -- The library's built-in overlay uses position:absolute with no positioned
  -- ancestor, so it covers the initial containing block (viewport at document
  -- origin) instead of the calendar area.  On scrollable Quarto pages the
  -- overlay therefore misses clicks near the calendar.  A document-level
  -- mousedown listener works around this reliably.
  table.insert(html_parts, [[
  document.addEventListener('mousedown', function(e) {
    var root = document.getElementById(]] .. utils.to_json(container_id) .. [[);
    if (!root) return;
    var popup = root.querySelector('.toastui-calendar-popup-container');
    if (!popup) return;
    if (popup.contains(e.target)) return;
    try { cal.getStoreDispatchers().popup.hideAllPopup(); } catch(ex) {}
  });]])

  if show_nav then
    table.insert(html_parts, nav_js(container_id))
  end

  table.insert(html_parts, '})();')
  table.insert(html_parts, '</script>')

  return pandoc.RawBlock("html", table.concat(html_parts, "\n"))
end

return M
