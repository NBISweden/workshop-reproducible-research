function Meta(meta)
  meta['quarto_version'] = tostring(quarto.version)
  meta['current_year'] = os.date("%Y")
  meta['current_date'] = os.date("%d-%m-%Y")
  meta['current_time'] = os.date("%H:%M:%S")
  meta['output-dir'] = quarto.project.output_directory

  local project_directory = quarto.project.directory or "."
  local quarto_config = io.open(project_directory .. "/_quarto.yml", "r")
  local blocks = {project = {}, website = {}, format = {}}
  local current_block = nil
  local stack = nil

  if quarto_config then
    for line in quarto_config:lines() do
      local block_name = line:match("^([%w%-]+):%s*$")

      if block_name and blocks[block_name] then
        current_block = block_name
        stack = {{indent = -1, value = blocks[current_block]}}
      elseif current_block and line:match("^%S") then
        current_block = nil
        stack = nil
      elseif current_block then
        local indentation = #(line:match("^(%s*)") or "")
        local key, value = line:match("^%s*([%w%-]+):%s*(.-)%s*$")

        if key then
          while stack[#stack].indent >= indentation do
            table.remove(stack)
          end

          local parent = stack[#stack].value
          if value == "" then
            parent[key] = {}
            table.insert(stack, {indent = indentation, value = parent[key]})
          elseif not value:match("^[{|]") then
            value = value:gsub('^[\"\']', ""):gsub('[\"\']$', "")
            parent[key] = value
          end
        end
      end
    end
    quarto_config:close()
  end

  meta.project = blocks.project
  meta.website = blocks.website
  meta.format = blocks.format

  return meta
end
