--- @module render
--- Team rendering helpers for HTML outputs.

local function load_module(path)
  return require(quarto.utils.resolve_path(path):gsub("%.lua$", ""))
end

local utils = load_module("_modules/utils.lua")

local M = {}

--- Emit an extension error in the current output format.
--- @param msg string
--- @param format string
--- @return pandoc.RawInline|pandoc.Strong
function M.error_inline(msg, format)
  if format == "html" then
    return pandoc.RawInline(
      "html",
      '<div class="team-error"><strong>Team Error:</strong> ' .. utils.escape_html_attr(msg) .. "</div>"
    )
  end

  return pandoc.Strong({ pandoc.Str("Team Error: " .. msg) })
end

--- Build an HTML error block for an invalid item.
--- @param msg string
--- @return pandoc.RawBlock
local function html_item_error(msg)
  return pandoc.RawBlock(
    "html",
    '<div class="team-error"><strong>Team Error:</strong> ' .. utils.escape_html_attr(msg) .. "</div>"
  )
end

--- Validate one team item.
--- @param item table
--- @param user_label string
--- @param index integer
--- @return string|nil
local function validate_item(item, user_label, index)
  local missing = {}
  if utils.is_empty(item.name) then table.insert(missing, "name") end
  if utils.is_empty(item.image) then table.insert(missing, "image") end
  if #missing == 0 then return nil end

  return string.format("'%s': Item %d is missing '%s'.", user_label, index, table.concat(missing, "' and '"))
end

--- Render team layout as HTML blocks.
--- @param team_dom_id string
--- @param team_items table
--- @param user_label string
--- @param name_class string|nil
--- @param image_class string|nil
--- @param description_class string|nil
--- @return pandoc.Blocks
function M.html(team_dom_id, team_items, user_label, name_class, image_class, description_class)
  local blocks = pandoc.Blocks({})

  local name_class_extra = (name_class ~= nil and name_class ~= "") and (" " .. utils.escape_html_attr(name_class)) or ""
  local image_class_attr = (image_class ~= nil and image_class ~= "") and (' class="' .. utils.escape_html_attr(image_class) .. '"') or ""
  local desc_class_extra = (description_class ~= nil and description_class ~= "") and (" " .. utils.escape_html_attr(description_class)) or ""

  blocks:insert(pandoc.RawBlock(
    "html",
    '<div id="' .. utils.escape_html_attr(team_dom_id) .. '" class="quarto-team" data-team-label="'
      .. utils.escape_html_attr(user_label) .. '">\n<div class="team-parent team-list">'
  ))

  for i, item in ipairs(team_items) do
    local item_error = validate_item(item, user_label, i)
    if item_error ~= nil then
      blocks:insert(html_item_error(item_error))
    else
      local item_dom_id = utils.build_item_id(team_dom_id, item, i)
      local name_text = pandoc.utils.stringify(item.name)
      local image_src = utils.escape_html_attr(pandoc.utils.stringify(item.image))
      local image_alt = utils.escape_html_attr(name_text)
      local image_url = item.image_url and pandoc.utils.stringify(item.image_url) or ""
      local name_url = item.name_url and pandoc.utils.stringify(item.name_url) or ""

      blocks:insert(pandoc.RawBlock(
        "html",
        '<article id="' .. utils.escape_html_attr(item_dom_id) .. '" class="team-child team-card">\n'
          .. '<div class="team-media">'
      ))

      if image_url ~= "" then
        blocks:insert(pandoc.RawBlock(
          "html",
          '<a class="team-image-link" href="' .. utils.escape_html_attr(image_url) .. '">'
        ))
      end

      blocks:insert(pandoc.RawBlock(
        "html",
        '<img src="' .. image_src .. '" alt="' .. image_alt .. '"' .. image_class_attr .. ' loading="lazy">'
      ))

      if image_url ~= "" then
        blocks:insert(pandoc.RawBlock("html", "</a>"))
      end

      blocks:insert(pandoc.RawBlock("html", "</div>"))

      blocks:insert(pandoc.RawBlock("html", '<div class="team-body">'))

      if name_url ~= "" then
        blocks:insert(pandoc.RawBlock(
          "html",
          '<a class="team-name-link" href="' .. utils.escape_html_attr(name_url) .. '">'
        ))
      end

      blocks:insert(pandoc.RawBlock("html", '<div class="name team-name' .. name_class_extra .. '">'))
      blocks:insert(pandoc.Plain(utils.content_to_inlines(item.name)))
      blocks:insert(pandoc.RawBlock("html", "</div>"))

      if name_url ~= "" then
        blocks:insert(pandoc.RawBlock("html", "</a>"))
      end

      if not utils.is_empty(item.description) then
        blocks:insert(pandoc.RawBlock("html", '<div class="description team-description' .. desc_class_extra .. '">'))
        blocks:extend(utils.content_to_blocks(item.description))
        blocks:insert(pandoc.RawBlock("html", "</div>"))
      end

      blocks:insert(pandoc.RawBlock("html", "</div>"))

      blocks:insert(pandoc.RawBlock("html", "</article>"))
    end
  end

  blocks:insert(pandoc.RawBlock("html", "</div>\n</div>"))
  return blocks
end

return M