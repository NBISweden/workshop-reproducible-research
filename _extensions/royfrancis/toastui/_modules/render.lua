--- Rendering utilities for calendar HTML and JavaScript payload.
--- @module toastui._modules.render

local utils = require('./utils')

local M = {}
local calendar_counter = 0

--- TOAST UI Calendar's own pristine default for `theme.common.backgroundColor`
--- (confirmed by inspecting a freshly-constructed instance with no custom
--- theme). Used as the fallback below so every calendar gets an explicit,
--- CSS-enforced value instead of silently reading whatever a *different*
--- calendar last wrote to it (see the comment on `theme_override_css`).
local DEFAULT_COMMON_BACKGROUND_COLOR = "white"

--- Scoped CSS reinforcement for `theme.common.backgroundColor`. TOAST UI
--- Calendar keeps this in state shared across every calendar instance on
--- the page: setting a custom background on one instance eventually
--- overwrites the inline style of every other instance too (including ones
--- built earlier), with whichever instance last touched the shared state
--- "winning" for the whole page. Emitting this rule — with the instance's
--- own requested color, or TOAST UI's own default when it didn't request
--- one — for *every* calendar, scoped by container id with `!important`,
--- beats that plain (non-`!important`) inline style regardless of what
--- other calendars on the page do or don't customize. Other theme fields
--- remain plain pass-through and stay subject to that upstream limitation
--- (see README's Limitations section).
--- @param container_id string
--- @param theme table|nil
--- @return string
local function theme_override_css(container_id, theme)
  local id = utils.escape_html_attr(container_id)
  local background_color = DEFAULT_COMMON_BACKGROUND_COLOR
  if type(theme) == "table" and type(theme.common) == "table" and theme.common.backgroundColor ~= nil then
    background_color = theme.common.backgroundColor
  end
  return '#' .. id .. ' .toastui-calendar-layout { background-color: ' .. utils.escape_html_attr(background_color) .. ' !important; }'
end
local ERROR_BOX_STYLE = "padding: 0.75rem 1rem; border: 1px solid #f5c6cb; border-radius: 6px; color: #721c24; background: #f8d7da;"
local NATIVE_DETAIL_ITEMS = {
  location = true,
  recurrenceRule = true,
  attendees = true,
  state = true,
  calendar = true,
  body = true,
}

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
    local escaped = utils.escape_html_attr(raw):gsub("\n", "<br>")
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
--- @param time_format string|nil
--- @param event_detail_items table|nil
--- @param popup_detail_items table|nil
--- @param auto_hour_range boolean
--- @return PandocRawBlock
function M.render_calendar_block(opts, calendars, events, initial_date, show_nav, height, timegrid_height, time_format, event_detail_items, popup_detail_items, auto_hour_range)
  local container_id = next_calendar_id()
  local html_parts = {}
  local has_custom_popup_details = false

  local tg = timegrid_height or "200%"
  table.insert(html_parts, '<style>#' .. utils.escape_html_attr(container_id) .. ' .toastui-calendar-timegrid { height: ' .. utils.escape_html_attr(tg) .. '; min-height: unset; }</style>')

  table.insert(html_parts, '<style>' .. theme_override_css(container_id, opts.theme) .. '</style>')

  if show_nav then
    local nav = nav_html(container_id)
    for _, line in ipairs(nav) do
      table.insert(html_parts, line)
    end
  else
    table.insert(html_parts, '<div class="toastui-calendar-wrapper">')
  end

  local detail_popup_classes = {}
  if opts.useDetailPopup then
    table.insert(detail_popup_classes, "toastui-calendar-detail-popup-enabled")

    local selected_detail_items = {}
    for _, item in ipairs(popup_detail_items or {}) do
      selected_detail_items[item] = true
    end
    local detail_item_classes = {
      { name = "location", class = "toastui-calendar-detail-exclude-location" },
      { name = "recurrenceRule", class = "toastui-calendar-detail-exclude-recurrence-rule" },
      { name = "attendees", class = "toastui-calendar-detail-exclude-attendees" },
      { name = "state", class = "toastui-calendar-detail-exclude-state" },
      { name = "calendar", class = "toastui-calendar-detail-exclude-calendar" },
      { name = "body", class = "toastui-calendar-detail-exclude-body" },
    }
    for _, item in ipairs(detail_item_classes) do
      if not selected_detail_items[item.name] then
        table.insert(detail_popup_classes, item.class)
      end
    end
    for _, item in ipairs(popup_detail_items or {}) do
      if not NATIVE_DETAIL_ITEMS[item] then
        has_custom_popup_details = true
        table.insert(detail_popup_classes, "toastui-calendar-detail-has-custom")
        break
      end
    end
  end
  table.insert(html_parts, '<div id="' .. utils.escape_html_attr(container_id) .. '" class="' .. table.concat(detail_popup_classes, " ") .. '" style="height: ' .. utils.escape_html_attr(height) .. ';"></div>')
  table.insert(html_parts, '</div>')

  table.insert(html_parts, '<script>')
  table.insert(html_parts, '(function() {')
  table.insert(html_parts, '  var opts = ' .. utils.to_json(opts) .. ';')

  if calendars then
    table.insert(html_parts, '  opts.calendars = ' .. utils.to_json(calendars) .. ';')
  end

  if time_format then
    table.insert(html_parts, '  var timeFormat = ' .. utils.to_json(time_format) .. ';')
    table.insert(html_parts, '  var eventDetailItems = ' .. utils.to_json(event_detail_items or {}) .. ';')
    table.insert(html_parts, '  var popupDetailItems = ' .. utils.to_json(popup_detail_items or {}) .. ';')
    table.insert(html_parts, [[  opts.template = opts.template || {};
  var detailItemIcons = {
    location: 'toastui-calendar-ic-location-b',
    recurrenceRule: 'toastui-calendar-ic-repeat-b',
    attendees: 'toastui-calendar-ic-user-b',
    state: 'toastui-calendar-ic-state-b',
    calendar: 'toastui-calendar-calendar-dot'
  };
  var nativeDetailItems = {
    location: true,
    recurrenceRule: true,
    attendees: true,
    state: true,
    calendar: true,
    body: true
  };
  function pad(value) {
    return String(value).padStart(2, '0');
  }
  function formatTime(value) {
    var hours = value.getHours();
    var minutes = pad(value.getMinutes());
    if (timeFormat === '24h') {
      return pad(hours) + ':' + minutes;
    }
    return (hours % 12 || 12) + ':' + minutes + (hours < 12 ? ' am' : ' pm');
  }
  function formatDate(value) {
    return value.getFullYear() + '.' + pad(value.getMonth() + 1) + '.' + pad(value.getDate());
  }
  function isSameDate(left, right) {
    return left.getFullYear() === right.getFullYear() &&
      left.getMonth() === right.getMonth() &&
      left.getDate() === right.getDate();
  }
  function escapeHtml(value) {
    var replacements = { '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' };
    return String(value == null ? '' : value).replace(/[&<>"']/g, function(character) {
      return replacements[character];
    });
  }
  function eventDetailValue(event, item) {
    if (item === 'calendar') {
      var calendars = opts.calendars || [];
      for (var index = 0; index < calendars.length; index++) {
        if (calendars[index].id === event.calendarId) return calendars[index].name || '';
      }
      return '';
    }

    var value = event[item];
    if (value == null && event.raw) value = event.raw[item];
    if (Array.isArray(value)) return value.join(', ');
    return value == null ? '' : value;
  }
  function renderEventDetails(event) {
    var details = eventDetailItems.map(function(item) {
      var value = eventDetailValue(event, item);
      if (value === '') return '';
      var iconClass = Object.prototype.hasOwnProperty.call(detailItemIcons, item) ? detailItemIcons[item] : '';
      var iconStyle = item === 'calendar' && event.backgroundColor ? ' style="background-color:' + escapeHtml(event.backgroundColor) + '"' : '';
      var icon = iconClass ? '<span class="toastui-calendar-icon toastui-calendar-event-detail-icon ' + iconClass + '"' + iconStyle + '></span>' : '';
      var label = Object.prototype.hasOwnProperty.call(nativeDetailItems, item) ? '' : '<strong>' + escapeHtml(item) + ':</strong> ';
      return '<span class="toastui-calendar-event-detail-item">' + icon + label + escapeHtml(value) + '</span>';
    }).join('');
    return details ? '<span class="toastui-calendar-event-detail-items">' + details + '</span>' : '';
  }
  function renderCustomPopupDetails(event) {
    var customItems = popupDetailItems.filter(function(item) {
      return !Object.prototype.hasOwnProperty.call(nativeDetailItems, item);
    });
    if (customItems.length === 0) return true;

    var section = container.querySelector('.toastui-calendar-section-detail');
    if (!section) return false;
    section.querySelectorAll('.toastui-calendar-custom-detail-item').forEach(function(element) {
      element.remove();
    });
    customItems.forEach(function(item) {
      var value = eventDetailValue(event, item);
      if (value === '') return;
      var row = document.createElement('div');
      row.className = 'toastui-calendar-detail-item toastui-calendar-custom-detail-item';
      var label = document.createElement('strong');
      label.textContent = item + ': ';
      row.appendChild(label);
      row.appendChild(document.createTextNode(String(value)));
      section.appendChild(row);
    });
    return true;
  }

  opts.template.time = function(event) {
    var title = escapeHtml(event.title);
    var titleAndTime = event.start ? '<strong>' + formatTime(event.start) + '</strong>&nbsp;' + title : title;
    return titleAndTime + renderEventDetails(event);
  };
  opts.template.timegridDisplayPrimaryTime = function(model) {
    return formatTime(model.time);
  };
  opts.template.timegridDisplayTime = function(model) {
    return formatTime(model.time);
  };
  opts.template.timegridNowIndicatorLabel = function(model) {
    return formatTime(model.time);
  };
  opts.template.popupDetailDate = function(event) {
    var sameDate = isSameDate(event.start, event.end);
    if (event.isAllday) {
      return formatDate(event.start) + (sameDate ? '' : ' - ' + formatDate(event.end));
    }

    var end = (sameDate ? '' : formatDate(event.end) + ' ') + formatTime(event.end);
    return formatDate(event.start) + ' ' + formatTime(event.start) + ' - ' + end;
  };
  Object.keys(opts.template).forEach(function(key) {
    if (typeof opts.template[key] !== 'string') return;
    var source = opts.template[key];
    opts.template[key] = function(model) {
      return source.replace(/\$\{\s*([\w.]+)\s*\}/g, function(_, path) {
        var value = path.split('.').reduce(function(acc, part) {
          return acc == null ? acc : acc[part];
        }, model);
        return value == null ? '' : escapeHtml(value);
      });
    };
  });]])
  end

  if events then
    -- Patch events so that detail-popup sections only appear for
    -- user-provided properties.  The library internally defaults
    -- state → "Busy" and attendees → [], both truthy, causing those
    -- sections to always render.  Setting explicit falsy values
    -- before createEvents() prevents those defaults from kicking in.
    --
    -- Built and normalized here, before the calendar is constructed,
    -- so autoHourRange (below) can measure already timezone-resolved
    -- event times and fold week.hourStart/hourEnd into opts before
    -- new tui.Calendar(opts) runs.
    table.insert(html_parts, '  var __ev = ' .. utils.to_json(events) .. ';')
    table.insert(html_parts, [[  var __sourceTz = (opts.timezone && opts.timezone.zones && opts.timezone.zones[0] && opts.timezone.zones[0].timezoneName) || null;
  var __naiveDateTime = /^(\d{4})-(\d{2})-(\d{2})T(\d{2}):(\d{2})(?::(\d{2}))?$/;
  function __tzOffsetMinutes(date, timeZone) {
    var parts = {};
    new Intl.DateTimeFormat('en-US', {
      timeZone: timeZone, hourCycle: 'h23',
      year: 'numeric', month: '2-digit', day: '2-digit',
      hour: '2-digit', minute: '2-digit', second: '2-digit'
    }).formatToParts(date).forEach(function(p) { parts[p.type] = p.value; });
    var asUTC = Date.UTC(parts.year, parts.month - 1, parts.day, parts.hour, parts.minute, parts.second);
    return (asUTC - date.getTime()) / 60000;
  }
  function __wallTimeToUtcIso(value, timeZone) {
    var match = __naiveDateTime.exec(value);
    if (!match) return null;
    var guess = Date.UTC(Number(match[1]), Number(match[2]) - 1, Number(match[3]), Number(match[4]), Number(match[5]), match[6] ? Number(match[6]) : 0);
    var offset = __tzOffsetMinutes(new Date(guess), timeZone);
    var utcMs = guess - offset * 60000;
    var refined = __tzOffsetMinutes(new Date(utcMs), timeZone);
    if (refined !== offset) utcMs = guess - refined * 60000;
    return new Date(utcMs).toISOString();
  }
  function __resolveEventTime(value) {
    if (__sourceTz && typeof value === 'string') {
      var converted = __wallTimeToUtcIso(value, __sourceTz);
      if (converted) return converted;
    }
    return value;
  }
  __ev.forEach(function(e, index) {
    if (e.id == null || e.id === '') e.id = ']] .. container_id .. [[-event-' + index;
    var source = Object.assign({}, e, e.raw && typeof e.raw === 'object' ? e.raw : {});
    e.raw = source;
    if (!e.state) e.state = '';
    if (!e.attendees || (Array.isArray(e.attendees) && e.attendees.length === 0)) e.attendees = null;
    if (!e.isAllday && e.category !== 'allday') {
      if (e.start != null) e.start = __resolveEventTime(e.start);
      if (e.end != null) e.end = __resolveEventTime(e.end);
    }
  });]])
    if auto_hour_range then
      table.insert(html_parts, [[  (function() {
    var weekOpts = opts.week || (opts.week = {});
    var hasHourStart = weekOpts.hourStart != null;
    var hasHourEnd = weekOpts.hourEnd != null;
    if (hasHourStart && hasHourEnd) return;
    var minHour = null, maxHour = null;
    // Shifts an instant by the delta between __sourceTz's offset and this
    // browser's own offset, so the native getters below read __sourceTz's
    // wall-clock digits (matching how TOAST UI actually draws the event)
    // instead of the viewer's local digits.
    function __zonedForDigits(d) {
      if (!__sourceTz) return d;
      return new Date(d.getTime() + (__tzOffsetMinutes(d, __sourceTz) + d.getTimezoneOffset()) * 60000);
    }
    __ev.forEach(function(e) {
      if (e.isAllday || e.category === 'allday') return;
      if (e.start == null || e.end == null) return;
      var start = new Date(e.start);
      var end = new Date(e.end);
      if (isNaN(start.getTime()) || isNaN(end.getTime())) return;
      start = __zonedForDigits(start);
      end = __zonedForDigits(end);
      var sameDay = start.getFullYear() === end.getFullYear() &&
        start.getMonth() === end.getMonth() &&
        start.getDate() === end.getDate();
      if (!sameDay) {
        minHour = 0;
        maxHour = 24;
        return;
      }
      var startHour = start.getHours() + start.getMinutes() / 60 + start.getSeconds() / 3600;
      var endHour = end.getHours() + end.getMinutes() / 60 + end.getSeconds() / 3600;
      if (minHour == null || startHour < minHour) minHour = startHour;
      if (maxHour == null || endHour > maxHour) maxHour = endHour;
    });
    if (minHour == null || maxHour == null) return;
    if (!hasHourStart) weekOpts.hourStart = minHour;
    if (!hasHourEnd) weekOpts.hourEnd = maxHour;
  })();]])
    end
  end

  table.insert(html_parts, '  var container = document.getElementById(' .. utils.to_json(container_id) .. ');')
  table.insert(html_parts, '  var cal = new tui.Calendar(container, opts);')
  if has_custom_popup_details then
    table.insert(html_parts, [[  cal.on('clickEvent', function(info) {
    var observer = new MutationObserver(function() {
      if (renderCustomPopupDetails(info.event)) observer.disconnect();
    });
    observer.observe(container, { childList: true, subtree: true });
    window.setTimeout(function() {
      renderCustomPopupDetails(info.event);
      observer.disconnect();
    }, 100);
  });]])
  end
  table.insert(html_parts, [[  function fitShortTimedEvent(event) {
    if (!event.start || !event.end || event.isAllday || event.category === 'allday') return;

    var duration = event.end.getTime() - event.start.getTime();
    var minimumDuration = 30 * 60 * 1000;
    if (duration <= 0 || duration >= minimumDuration) return;

    var weekOptions = opts.week || {};
    var hourStart = weekOptions.hourStart == null ? 0 : Number(weekOptions.hourStart);
    var hourEnd = weekOptions.hourEnd == null ? 24 : Number(weekOptions.hourEnd);
    var visibleDuration = (hourEnd - hourStart) * 60 * 60 * 1000;
    if (visibleDuration <= 0) return;

    var height = duration / visibleDuration * 100;
    if (event.id == null || event.id === '') return;
    var eventId = String(event.id);
    container.querySelectorAll('.toastui-calendar-event-time').forEach(function(element) {
      if (element.getAttribute('data-event-id') !== eventId) return;
      element.style.height = 'max(1px, calc(' + height + '% - 2px))';
      element.style.minHeight = '0';
      element.style.overflow = 'hidden';
      var content = element.querySelector('.toastui-calendar-event-time-content');
      if (content) {
        content.style.minHeight = '0';
        content.style.overflow = 'hidden';
      }
    });
  }
  cal.on('afterRenderEvent', fitShortTimedEvent);]])

  if events then
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
