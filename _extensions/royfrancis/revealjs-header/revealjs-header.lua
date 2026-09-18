-- quarto revealjs-header extension

--- Register JavaScript and CSS assets for this extension.
--- This is safe to call once for revealjs documents.
local function ensure_html_dependencies()
  quarto.doc.add_html_dependency({
    name = "revealjs-header",
    scripts = {
      { path = "revealjs-header.js", attribs = { defer = "true" } }
    },
    stylesheets = { "revealjs-header.css" }
  })
end

--- Convert a metadata value to a non-empty string.
--- @param value any Metadata value.
--- @return string|nil
local function to_non_empty_string(value)
  if value == nil then
    return nil
  end

  local as_string = pandoc.utils.stringify(value)
  if as_string == nil then
    return nil
  end

  if as_string:match("^%s*$") then
    return nil
  end

  return as_string
end

--- Create a logo image inline.
--- @param path string Image path.
--- @param height string|nil CSS height (for example "50px" or "2em").
--- @return Inline
local function make_image(path, height)
  local image_attr = pandoc.Attr("", {}, {})

  if height ~= nil then
    image_attr = pandoc.Attr("", {}, {
      { "style", "height:" .. height .. ";max-width:none;" }
    })
  end

  return pandoc.Image("", path, "", image_attr)
end

--- Wrap an inline in a link when URL exists.
--- @param inline Inline The inline to wrap.
--- @param url string|nil URL target.
--- @return Inline
local function with_optional_link(inline, url)
  if url == nil then
    return inline
  end

  return pandoc.Link(inline, url)
end

--- Build one logo container as a Div block.
--- @param path string|nil Logo image path.
--- @param url string|nil Optional link for the logo.
--- @param height string|nil Optional CSS height.
--- @param class_name string CSS class to apply to the container.
--- @return Block
local function make_logo_div(path, url, height, class_name)
  if path == nil then
    -- Keep an empty placeholder so the two-column layout remains stable.
    return pandoc.Div({}, { class = class_name })
  end

  local image = make_image(path, height)
  local content = with_optional_link(image, url)
  return pandoc.Div({ pandoc.Plain({ content }) }, { class = class_name })
end

--- Read all extension options from document metadata.
--- @param meta Meta
--- @return table
local function read_options(meta)
  return {
    left_path = to_non_empty_string(meta["header-logo-left"]),
    left_url = to_non_empty_string(meta["header-logo-left-url"]),
    left_height = to_non_empty_string(meta["header-logo-left-height"]),
    right_path = to_non_empty_string(meta["header-logo-right"]),
    right_url = to_non_empty_string(meta["header-logo-right-url"]),
    right_height = to_non_empty_string(meta["header-logo-right-height"])
  }
end

if quarto.doc.is_format("revealjs") then
  ensure_html_dependencies()
end

--- Pandoc document filter entrypoint.
--- @param doc Pandoc
--- @return Pandoc
function Pandoc(doc)
  if not quarto.doc.is_format("revealjs") then
    return doc
  end

  local options = read_options(doc.meta)
  local left_logo = make_logo_div(
    options.left_path,
    options.left_url,
    options.left_height,
    "header-logo-left"
  )
  local right_logo = make_logo_div(
    options.right_path,
    options.right_url,
    options.right_height,
    "header-logo-right"
  )

  local header = pandoc.Div(
    { left_logo, right_logo },
    { class = "reveal-header" }
  )

  table.insert(doc.blocks, header)
  return doc
end
