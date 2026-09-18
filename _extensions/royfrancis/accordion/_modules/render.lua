--- @module render
--- Accordion rendering helpers for html, pdf, typst, and fallback formats.

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
      '<div class="accordion-error"><strong>Accordion Error:</strong> ' .. msg .. "</div>"
    )
  end

  if format == "pdf" then
    return pandoc.RawInline(
      "latex",
      "{\\color{accordionerror}\\textbf{Accordion Error:} " .. utils.escape_latex(msg) .. "}"
    )
  end

  if format == "typst" then
    return pandoc.RawInline(
      "typst",
      '#text(fill: rgb("#721c24"))[*Accordion Error:* ' .. utils.escape_typst(msg) .. "]"
    )
  end

  return pandoc.Strong({ pandoc.Str("Accordion Error: " .. msg) })
end

--- Build an HTML error block for an invalid item.
--- @param msg string
--- @return pandoc.RawBlock
local function html_item_error(msg)
  return pandoc.RawBlock(
    "html",
    '<div class="accordion-error"><strong>Accordion Error:</strong> ' .. msg .. "</div>"
  )
end

--- Render accordion as Bootstrap HTML structure.
--- @param accordion_id string
--- @param accordion_items table
--- @param user_label string
--- @return pandoc.Blocks
function M.html(accordion_id, accordion_items, user_label)
  local blocks = pandoc.Blocks({})

  blocks:insert(pandoc.RawBlock("html", '<div id="' .. accordion_id .. '" class="accordion quarto-accordion">'))

  for i, item in ipairs(accordion_items) do
    local header = pandoc.utils.stringify(item.header or "")
    local body = pandoc.utils.stringify(item.body or "")

    local item_error = utils.item_validation_error(user_label, i, header, body)
    if item_error ~= nil then
      blocks:insert(html_item_error(item_error))
    else
      local collapse_id = utils.build_collapse_id(accordion_id, item, i)
      local collapsed = item.collapsed
      if collapsed == nil then collapsed = true end

      local collapse_class = collapsed and "collapsed" or ""
      local aria_expanded = collapsed and "false" or "true"
      local collapse_show = collapsed and "" or " show"

      blocks:insert(pandoc.RawBlock("html",
        '<div class="accordion-item">\n'
          .. '<div class="accordion-header" id="heading' .. collapse_id .. '">\n'
          .. '<button class="accordion-button ' .. collapse_class .. '" type="button" '
          .. 'data-bs-toggle="collapse" data-bs-target="#collapse' .. collapse_id .. '" '
          .. 'aria-expanded="' .. aria_expanded .. '" aria-controls="collapse' .. collapse_id .. '">\n'
          .. '<div class="accordion-header-content">'))

      blocks:extend(utils.content_to_blocks(item.header))

      blocks:insert(pandoc.RawBlock("html",
        '</div>\n</button>\n</div>\n'
          .. '<div id="collapse' .. collapse_id .. '" class="accordion-collapse collapse' .. collapse_show .. '" '
          .. 'aria-labelledby="heading' .. collapse_id .. '" data-bs-parent="#' .. accordion_id .. '">\n'
          .. '<div class="accordion-body">\n'
          .. '<div class="accordion-body-content">'))

      blocks:extend(utils.content_to_blocks(item.body))

      blocks:insert(pandoc.RawBlock("html", "</div>\n</div>\n</div>\n</div>"))
    end
  end

  blocks:insert(pandoc.RawBlock("html", "</div>"))
  return blocks
end

--- Render accordion as card-style LaTeX blocks for PDF output.
--- @param accordion_items table
--- @param user_label string
--- @return pandoc.Blocks
function M.pdf(accordion_items, user_label)
  local blocks = pandoc.Blocks({})

  for i, item in ipairs(accordion_items) do
    local header = pandoc.utils.stringify(item.header or "")
    local body = pandoc.utils.stringify(item.body or "")

    local item_error = utils.item_validation_error(user_label, i, header, body)
    if item_error ~= nil then
      blocks:insert(pandoc.RawBlock(
        "latex",
        "{\\color{accordionerror}\\textbf{Accordion Error:} " .. utils.escape_latex(item_error) .. "}"
      ))
    else
      blocks:insert(pandoc.RawBlock(
        "latex",
        "\\begin{tcolorbox}[colback=white, colframe=black!20, boxrule=0.5pt, arc=3pt]"
      ))
      blocks:extend(utils.make_bold(utils.content_to_blocks(item.header)))
      blocks:insert(pandoc.RawBlock("latex", "\\tcblower"))
      blocks:extend(utils.content_to_blocks(item.body))
      blocks:insert(pandoc.RawBlock("latex", "\\end{tcolorbox}"))
    end
  end

  return blocks
end

--- Render accordion as card-style Typst blocks.
--- @param accordion_items table
--- @param user_label string
--- @return pandoc.Blocks
function M.typst(accordion_items, user_label)
  local blocks = pandoc.Blocks({})

  for i, item in ipairs(accordion_items) do
    local header = pandoc.utils.stringify(item.header or "")
    local body = pandoc.utils.stringify(item.body or "")

    local item_error = utils.item_validation_error(user_label, i, header, body)
    if item_error ~= nil then
      blocks:insert(pandoc.RawBlock(
        "typst",
        '#text(fill: rgb("#721c24"))[*Accordion Error:* ' .. utils.escape_typst(item_error) .. "]"
      ))
    else
      blocks:insert(pandoc.RawBlock(
        "typst",
        "#block(width: 100%, stroke: 0.5pt + luma(180), radius: 4pt, inset: 0pt, below: 8pt, breakable: false)[\n"
          .. "#block(inset: 10pt, width: 100%, below: 0pt, above: 0pt)["
      ))
      blocks:extend(utils.make_bold(utils.content_to_blocks(item.header)))
      blocks:insert(pandoc.RawBlock(
        "typst",
        "]\n#line(length: 100%, stroke: 0.3pt + luma(200))\n"
          .. "#block(inset: 10pt, width: 100%, below: 0pt, above: 0pt)["
      ))
      blocks:extend(utils.content_to_blocks(item.body))
      blocks:insert(pandoc.RawBlock("typst", "]\n]"))
    end
  end

  return blocks
end

--- Render accordion in plain Pandoc blocks as a fallback.
--- @param accordion_items table
--- @param user_label string
--- @return pandoc.Blocks
function M.fallback(accordion_items, user_label)
  local blocks = pandoc.Blocks({})

  for i, item in ipairs(accordion_items) do
    local header = pandoc.utils.stringify(item.header or "")
    local body = pandoc.utils.stringify(item.body or "")

    local item_error = utils.item_validation_error(user_label, i, header, body)
    if item_error ~= nil then
      blocks:insert(pandoc.Para({ pandoc.Strong({ pandoc.Str("Accordion Error: " .. item_error) }) }))
    else
      blocks:extend(utils.make_bold(utils.content_to_blocks(item.header)))
      blocks:extend(utils.content_to_blocks(item.body))
    end

    if i < #accordion_items then
      blocks:insert(pandoc.HorizontalRule())
    end
  end

  return blocks
end

return M