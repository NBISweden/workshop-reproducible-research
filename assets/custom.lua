local function meta_inlines(value)
  return pandoc.MetaInlines({pandoc.Str(tostring(value))})
end

local UNSUPPORTED_VALUE = "[unsupported yaml]"

function Meta(meta)
  meta['quarto_version'] = meta_inlines(quarto.version)
  meta['current_year'] = meta_inlines(os.date("%Y"))
  meta['current_date'] = meta_inlines(os.date("%d-%m-%Y"))
  meta['current_time'] = meta_inlines(os.date("%H:%M:%S"))
  meta['output-dir-path'] = meta_inlines(quarto.project.output_directory)
  meta['output-dir'] = meta_inlines(quarto.project.output_directory:match("([^/\\]+)[/\\]*$"))

  local project_directory = quarto.project.directory or "."
  local quarto_config = io.open(project_directory .. "/_quarto.yml", "r")
  local blocks = {project = {}, website = {}, format = {}}
  local current_block = nil
  local stack = nil
  local skip_until = nil

  if quarto_config then
    for line in quarto_config:lines() do
      local content = line:match("^%s*(.-)%s*$")

      if content ~= "" and content:sub(1, 1) ~= "#" then
        local block_name = line:match("^([%w%-]+):%s*$")

        if block_name and blocks[block_name] then
          current_block = block_name
          stack = {{indent = -1, value = blocks[current_block]}}
          skip_until = nil
        elseif current_block and line:match("^%S") then
          current_block = nil
          stack = nil
          skip_until = nil
        elseif current_block then
          local indentation = #(line:match("^(%s*)") or "")

          if not (skip_until ~= nil and indentation > skip_until) then
            skip_until = nil

            if content:sub(1, 1) == "-" then
              -- block-style YAML list: unsupported, mark parent and skip its items
              skip_until = indentation
              local top = stack[#stack]
              if top.parent then
                top.parent[top.key] = meta_inlines(UNSUPPORTED_VALUE)
              end
            else
              local key, value = line:match("^%s*([%w%-]+):%s*(.-)%s*$")

              if key then
                while stack[#stack].indent >= indentation do
                  table.remove(stack)
                end

                local parent = stack[#stack].value
                if value == "" then
                  parent[key] = {}
                  table.insert(stack, {indent = indentation, value = parent[key], key = key, parent = parent})
                elseif value:match("^[%[{|>]") then
                  -- flow-style list/map or block scalar: unsupported
                  parent[key] = meta_inlines(UNSUPPORTED_VALUE)
                else
                  value = value:gsub('^[\"\']', ""):gsub('[\"\']$', "")
                  parent[key] = meta_inlines(value)
                end
              end
            end
          end
        end
      end
    end
    quarto_config:close()
  end

  meta.custom = blocks

  return meta
end
