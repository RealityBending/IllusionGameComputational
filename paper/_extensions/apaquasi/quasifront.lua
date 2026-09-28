-- apaquasi: the first page of a preprint, laid out the way a high-impact
-- journal lays out the first page of an article. The model is Nature Human
-- Behaviour.
--
-- Called by frontmatter.lua in journal mode for the pdf, in place of the
-- front matter apaquarto builds there. Everything else is left to apaquarto:
-- the blocks returned here are handed to formatlatex.lua under the same div
-- classes apaquarto's journal front matter uses, so the masthead still goes to
-- \twocolumn and the note still goes to the foot of the first column.
--
--   JournalMasthead  the band, the title, the sidebar, the byline (with the
--                    orcids), the affiliations and the abstract, spanning
--                    the page
--   JournalNote      the author note (the disclosures, the CRediT roles) and
--                    the correspondence address, at the foot of the first
--                    column
--
-- The fields it reads, all optional, sit under titlepage: in the yaml:
--
--   titlepage:
--     type: Research Article   # the wordmark under the bar (the default):
--                              # Review, Tutorial, Registered Report...
--     label: Preprint          # over the title, at the left (the default)
--     server: PsyArXiv         # the sidebar, and the foot of the page
--     doi: 10.31234/osf.io/abcde   # over the title, at the right
--     version: 1
--     status: Not peer reviewed
--     license: CC BY 4.0
--     color: 1B5E9F            # the accent, as hex
--     date-label: Posted       # what the date is called in the sidebar
--     badges:                  # open science badges: a url, or true
--       open-data: https://osf.io/abcde     # or open-data-pa
--       open-materials: https://github.com/...
--       open-code: https://github.com/...
--       preregistered: https://osf.io/xyz   # or preregistered-plus, and
--                                           # either with -tc, -de or -de-tc
--
-- and the document's own date, keywords and supplemental-materials. The
-- materials may be one url or a list of them; each is listed with the icon of
-- its repository as its bullet.
--
-- Latex is written in long strings, [[...]], which take a backslash as it is.

local M = {}

local List = require 'pandoc.List'
local utilsapa = require("utilsapa")
local stringify = utilsapa.stringify

local function raw(text) return pandoc.RawBlock("latex", text) end
local function rawi(text) return pandoc.RawInline("latex", text) end

-- A metadata value as inlines, or nil when it is absent or empty.
local function inlines(value)
  if value == nil then return nil end
  local kind = pandoc.utils.type(value)
  if kind == "Inlines" then
    if #value == 0 then return nil end
    return value
  end
  if kind == "Blocks" then
    return pandoc.utils.blocks_to_inlines(value)
  end
  local text = stringify(value)
  if text == "" or text == "false" then return nil end
  return pandoc.Inlines(text)
end

local function field(meta, name)
  local pp = meta.titlepage
  if pp == nil or pandoc.utils.type(pp) ~= "table" then return nil end
  return pp[name]
end

local function is_true(value)
  return value ~= nil and stringify(value) == "true"
end

local function join(parts, sep)
  local out = pandoc.Inlines({})
  for i, part in ipairs(parts) do
    if i > 1 then out:extend(pandoc.Inlines(sep)) end
    out:extend(part)
  end
  return out
end

-- Inlines wrapped in a latex command: \name{...}.
local function command(name, content)
  local out = pandoc.Inlines({ rawi([[\]] .. name .. "{") })
  out:extend(content)
  out:insert(rawi("}"))
  return out
end

-- The paragraphs of an abstract, from either of the shapes it arrives in: a
-- line block to a paragraph when written under a yaml block scalar, ordinary
-- paragraphs when written as a section of the document.
local function paragraphs(value)
  local out = List:new {}
  local kind = pandoc.utils.type(value)
  if kind == "Inlines" then
    out:insert(pandoc.Para(value))
  elseif kind == "Blocks" then
    for _, block in ipairs(value) do
      if block.t == "LineBlock" then
        for _, line in ipairs(block.content) do
          out:insert(pandoc.Para(line))
        end
      elseif block.t == "Para" or block.t == "Plain" then
        out:insert(pandoc.Para(block.content))
      elseif block.t ~= "Header" then
        out:insert(block)
      end
    end
  end
  return out
end

local function masked(meta)
  return is_true(meta.mask)
end

local function corresponding_author(meta)
  for _, a in ipairs(meta["by-author"] or {}) do
    if a.attributes and is_true(a.attributes.corresponding) then
      return a
    end
  end
end

-- The byline: every author, a superscript for each of their affiliations,
-- the orcid mark linked to their record, and an envelope after the
-- corresponding author.
local function byline(meta)
  local authors = meta["by-author"]
  if authors == nil or #authors == 0 or meta["suppress-author"] then return nil end
  local out = pandoc.Inlines({})
  for i, a in ipairs(authors) do
    if i > 1 then
      out:extend(pandoc.Inlines(i == #authors and " & " or ", "))
    end
    out:extend(inlines(a.apaauthordisplay) or pandoc.Inlines(""))
    local numbers = {}
    for _, aff in ipairs(a.affiliations or {}) do
      if aff.number then table.insert(numbers, stringify(aff.number)) end
    end
    local marks = pandoc.Inlines({})
    if #numbers > 0 and not meta["suppress-affiliation"] then
      marks:insert(pandoc.Superscript(table.concat(numbers, ",")))
    end
    if a.orcid and stringify(a.orcid) ~= "" and not meta["suppress-orcid"] then
      marks:insert(rawi([[\,\orcidlink{]] .. stringify(a.orcid) .. "}"))
    end
    if a.attributes and is_true(a.attributes.corresponding) then
      marks:insert(rawi([[\,\quasienvelope]]))
    end
    out:extend(marks)
  end
  return out
end

-- One affiliation, run in: its number, then department, institution, city
-- and country, the way the journals list them.
local function affiliation_line(aff)
  local parts = {}
  for _, key in ipairs({ "group", "department", "name", "address", "city", "region", "country" }) do
    local value = inlines(aff[key])
    if value then table.insert(parts, value) end
  end
  if #parts == 0 then return nil end
  local out = pandoc.Inlines({})
  if aff.number then out:insert(pandoc.Superscript(stringify(aff.number))) end
  out:extend(join(parts, ", "))
  out:extend(pandoc.Inlines("."))
  return out
end

local function is_picture(inline)
  if inline.t == "Image" then return inline.identifier ~= "orcid" end
  return inline.t == "Link" and #inline.content == 1 and inline.content[1].t == "Image"
end

-- A picture in a paragraph of the note (a logo in the gratitude note), which
-- apaquarto runs in with the rest of the disclosures, taken out to a line of
-- its own after the paragraph. In the narrow foot, run in, it split a line of
-- text in two.
local function pictures_apart(block)
  if block.t ~= "Para" then return { block } end
  local text, pictures = pandoc.Inlines({}), pandoc.Inlines({})
  for _, inline in ipairs(block.content) do
    if is_picture(inline) then pictures:insert(inline) else text:insert(inline) end
  end
  if #pictures == 0 then return { block } end
  while #text > 0 and (text[#text].t == "Space" or text[#text].t == "SoftBreak" or
      text[#text].t == "LineBreak") do
    text:remove()
  end
  local out = {}
  if #text > 0 then table.insert(out, pandoc.Para(text)) end
  table.insert(out, pandoc.Para(pictures))
  return out
end

-- The affiliations, one to a line, for under the byline.
local function affiliations(meta)
  if not meta.affiliations or meta["suppress-affiliation"] then return nil end
  local out = pandoc.Inlines({})
  for _, aff in ipairs(meta.affiliations) do
    local line = affiliation_line(aff)
    if line then
      if #out > 0 then out:insert(pandoc.LineBreak()) end
      out:extend(line)
    end
  end
  if #out == 0 then return nil end
  return out
end

-- An orcid line of apaquarto's author note: a name, the mark and the link.
-- The byline carries the mark, so the line is left out of the note.
local function is_orcid_line(block)
  local found = false
  if block.t == "Para" then
    block:walk {
      RawInline = function(r)
        if r.text:match("^\\orcidlink") then found = true end
      end,
      Image = function(img)
        if img.identifier == "orcid" then found = true end
      end
    }
  end
  return found
end

-- The author note as apaquarto builds it, less the orcid lines and the
-- address, which is written here in the journals' own way.
local function author_note_blocks(notes, corresponding)
  local kept = List:new {}
  for _, block in ipairs(notes or {}) do
    if not (corresponding and block == corresponding) and not is_orcid_line(block) then
      kept:extend(pictures_apart(block))
    end
  end
  return kept
end

-- The foot of the first column: the author note, and the address to write to.
local function footnote_blocks(meta, notes, corresponding)
  local out = List:new {}
  if masked(meta) then return out end

  out:extend(author_note_blocks(notes, corresponding))

  local note = meta["author-note"]
  local address
  if note and note["correspondence-note"] then
    address = inlines(note["correspondence-note"])
  else
    local a = corresponding_author(meta)
    if a and a.email and not meta["suppress-corresponding-email"] then
      local email = stringify(a.email)
      address = pandoc.Inlines({ rawi([[\quasienvelope\,]]) })
      address:extend(pandoc.Inlines("e-mail: "))
      address:insert(pandoc.Link(email, "mailto:" .. email))
    end
  end
  if address and not meta["suppress-corresponding-paragraph"] then
    out:insert(pandoc.Para(address))
  end
  return out
end

-- A link shown without its scheme, which is how the journals print one. In a
-- narrow column it may break after any slash.
local function short_link(url)
  local shown = url:gsub("^https?://", ""):gsub("^www%.", ""):gsub("/+$", "")
  local text = pandoc.Inlines({})
  local first = true
  for part in (shown .. "/"):gmatch("(.-)/") do
    if not first then
      text:insert(pandoc.Str("/"))
      text:insert(rawi([[\allowbreak{}]]))
    end
    if part ~= "" then text:insert(pandoc.Str(part)) end
    first = false
  end
  return pandoc.Inlines({ pandoc.Link(text, url) })
end

local function is_url(text)
  return text:match("^https?://%S+$") ~= nil
end

-- The icon a materials link is listed with: the repository's own where there
-- is one, a link otherwise.
local repository_icons = {
  { "github%.com", [[\faGithub]] },
  { "gitlab%.",    [[\faGitlab]] },
  { "osf%.io",     [[\quasiai{osf}]] },
  { "zenodo%.org", [[\quasiai{zenodo}]] },
  { "psyarxiv",    [[\quasiai{psyarxiv}]] },
  { "figshare",    [[\quasiai{figshare}]] },
}

local function repository_icon(url)
  local host = url:match("^https?://([^/]+)") or ""
  for _, pair in ipairs(repository_icons) do
    if host:match(pair[1]) then return pair[2] end
  end
  return [[\faLink]]
end

-- The materials, as a list with the repository's icon as each bullet. They
-- may be written as one value or a list; anything that is not a url is kept
-- as written, with the link icon.
local function materials_list(value)
  if value == nil or stringify(value) == "" then return nil end
  local items = {}
  if pandoc.utils.type(value) == "List" then
    for _, item in ipairs(value) do table.insert(items, item) end
  else
    table.insert(items, value)
  end
  local out = pandoc.Inlines({})
  for _, item in ipairs(items) do
    local text = stringify(item)
    local icon, shown = [[\faLink]], inlines(item)
    if is_url(text) then
      icon, shown = repository_icon(text), short_link(text)
    end
    if shown then
      out:insert(rawi([[\quasimaterial{]] .. icon .. "}{"))
      out:extend(shown)
      out:insert(rawi("}"))
    end
  end
  if #out == 0 then return nil end
  return out
end

-- The open science badges, in the order the Center for Open Science gives
-- them. Each is asked for by its key, with the url of what it certifies or
-- with true, and drawn from the badge files the Center publishes
-- (https://osf.io/tvyxz/), which ship in badges/: vector cut from its sheet
-- where there is one, its own png for the two the sheet leaves out. TC is
-- transparent changes, DE data exist.
local badge_kinds = {
  { "open-data",                "open-data.pdf" },
  { "open-data-pa",             "open-data-pa.pdf" },
  { "open-materials",           "open-materials.pdf" },
  { "open-code",                "open-code.png" },
  { "preregistered",            "preregistered.pdf" },
  { "preregistered-tc",         "preregistered-tc.pdf" },
  { "preregistered-de",         "preregistered-de.pdf" },
  { "preregistered-de-tc",      "preregistered-de-tc.pdf" },
  { "preregistered-plus",       "preregistered-plus.png" },
  { "preregistered-plus-tc",    "preregistered-plus-tc.pdf" },
  { "preregistered-plus-de",    "preregistered-plus-de.pdf" },
  { "preregistered-plus-de-tc", "preregistered-plus-de-tc.pdf" },
}

-- How many badges a row of the sidebar holds: three at the size they are set
-- at, and four, a little smaller, when there are more than three, so that a
-- fourth does not stand on a row of its own.
local function badges_per_row(count)
  if count <= 3 then return 3 end
  return 4
end

-- The width of one badge, as a fraction of the sidebar, for that many to a
-- row with the gaps between them (\quasibadgegap, 0.035 of the sidebar).
-- Rounded down, so that rounding never pushes the last of a row onto the next.
local function badge_width(count)
  local per_row = badges_per_row(count)
  local width = (1 - 0.035 * (per_row - 1)) / per_row
  return string.format("%.4f", math.floor(width * 10000) / 10000)
end

local function badges_row(meta)
  local asked = field(meta, "badges")
  if asked == nil or pandoc.utils.type(asked) ~= "table" then return nil end
  local shown = {}
  for _, kind in ipairs(badge_kinds) do
    local value = asked[kind[1]]
    local text = value ~= nil and stringify(value) or ""
    if text ~= "" and text ~= "false" then table.insert(shown, { kind, text }) end
  end
  if #shown == 0 then return nil end

  local out = pandoc.Inlines({ rawi([[\quasibadges{]] .. badge_width(#shown) .. "}") })
  for i, item in ipairs(shown) do
    local kind, text = item[1], item[2]
    if i > 1 then out:insert(rawi([[\quasibadgegap]])) end
    local file = "badges/" .. kind[2]
    file = utilsapa.extension_file_relative(file) or file
    local badge = rawi([[\quasibadge{]] .. file .. "}")
    if is_url(text) then
      out:insert(pandoc.Link({ badge }, text))
    else
      out:insert(badge)
    end
  end
  return out
end

-- The sidebar at the left of the byline and the abstract: what the journals
-- print there (the dates) and what a preprint has in their place. Each entry
-- is ruled off from the next.
local function sidebar_blocks(meta)
  local out = List:new {}
  local function entry(label, value)
    if value == nil then return end
    local content = pandoc.Inlines({})
    if label then content:extend(command("quasisidelabel", pandoc.Inlines(label))) end
    content:extend(value)
    out:insert(pandoc.Para(content))
    out:insert(raw([[\quasisiderule]]))
  end

  entry("Server", inlines(field(meta, "server")))
  entry(stringify(field(meta, "date-label"), "Posted"), inlines(meta.date))
  entry("Version", inlines(field(meta, "version")))
  entry(nil, inlines(field(meta, "status")))
  entry("License", inlines(field(meta, "license")))

  if meta.keywords and not meta["suppress-keywords"] then
    local words = {}
    if pandoc.utils.type(meta.keywords) == "List" then
      for _, k in ipairs(meta.keywords) do table.insert(words, inlines(k)) end
    else
      table.insert(words, inlines(meta.keywords))
    end
    if #words > 0 then entry("Keywords", join(words, ", ")) end
  end

  if not meta["suppress-supplemental-materials"] then
    entry("Materials", materials_list(meta["supplemental-materials"]))
  end

  local badges = badges_row(meta)
  if badges then
    out:insert(pandoc.Para(badges))
    out:insert(raw([[\quasisiderule]]))
  end
  return out
end

local function abstract_blocks(meta)
  local abstract = List:new {}
  if meta.apaabstract and #meta.apaabstract > 0 and not meta["suppress-abstract"] then
    abstract:extend(paragraphs(meta.apaabstract))
  end
  if meta["impact-statement"] and #meta["impact-statement"] > 0 and
      not meta["suppress-impact-statement"] then
    local heading = inlines(meta.language and meta.language["title-impact-statement"]) or
      pandoc.Inlines("Impact Statement")
    local statement = paragraphs(meta["impact-statement"])
    if #statement > 0 and statement[1].t == "Para" then
      local first = pandoc.Inlines({ rawi([[{\sffamily\bfseries\small ]]) })
      first:extend(heading)
      first:insert(rawi(".}"))
      first:insert(pandoc.Space())
      first:extend(statement[1].content)
      statement[1] = pandoc.Para(first)
    end
    abstract:extend(statement)
  end
  return abstract
end

-- The band, the title, then the sidebar beside the byline and the abstract.
local function masthead_blocks(meta)
  local out = List:new {}
  local server = inlines(field(meta, "server"))
  local date = inlines(meta.date)
  -- The wordmark says what kind of article this is, as a journal's section
  -- would; the label over the title, what stage it is at. Either is hidden
  -- by writing false.
  local kind = inlines(field(meta, "type"))
  if field(meta, "type") == nil then kind = pandoc.Inlines("Research Article") end
  local label = inlines(field(meta, "label"))
  if field(meta, "label") == nil then label = pandoc.Inlines("Preprint") end

  -- The accent and the foot of the first page are set from inside the box
  -- formatlatex.lua builds the masthead in, so both are made global.
  local color = field(meta, "color")
  if color then
    local hex = stringify(color):gsub("^#", "")
    out:insert(raw([[\xglobal\definecolor{quasiaccent}{HTML}{]] .. hex .. "}"))
  end
  local foot = List:new {}
  if label then foot:insert(label) end
  if server then foot:insert(server) end
  if date then foot:insert(date) end
  out:insert(pandoc.Plain(command("setquasifootline", join(foot, " | "))))

  -- The bar, and the wordmark under it.
  out:insert(raw([[\quasibar]]))
  if kind then
    out:insert(pandoc.Plain(command("quasiwordmark", kind)))
  end

  -- The label, and the doi opposite it.
  local doi = inlines(field(meta, "doi"))
  local typeline = command("quasitypeline", label or pandoc.Inlines({}))
  typeline:insert(rawi("{"))
  if doi then typeline:extend(short_link("https://doi.org/" .. stringify(doi))) end
  typeline:insert(rawi("}"))
  out:insert(pandoc.Plain(typeline))

  if meta.apatitledisplay and not meta["suppress-title"] then
    out:insert(raw([[\begin{quasititle}]]))
    out:insert(pandoc.Plain(meta.apatitledisplay))
    out:insert(raw([[\end{quasititle}]]))
  end

  local side = sidebar_blocks(meta)
  local authors = not masked(meta) and byline(meta)
  local places = not masked(meta) and affiliations(meta)
  local abstract = abstract_blocks(meta)
  if #side > 0 or authors or #abstract > 0 then
    out:insert(raw([[\quasisidebar]]))
    out:extend(side)
    out:insert(raw([[\quasimain]]))
    if authors then out:insert(pandoc.Plain(authors)) end
    if places then
      out:insert(raw([[\quasiaffiliations]]))
      out:insert(pandoc.Plain(places))
    end
    if #abstract > 0 then
      out:insert(raw([[\quasiabstract]]))
      out:extend(abstract)
    end
    out:insert(raw([[\quasimainend]]))
  end
  return out
end

-- APA's limit on the abstract. The first page has room for about this much
-- beside a sidebar of the usual length; apaquasi.tex scales down a front
-- matter taller than the page, and this says why it may have been.
local kAbstractWords = 250

local function warn_long_abstract(meta)
  if not meta.apaabstract or meta["suppress-abstract"] then return end
  local words = 0
  for _ in stringify(meta.apaabstract):gmatch("%S+") do words = words + 1 end
  if words > kAbstractWords then
    quarto.log.warning("The abstract is " .. words .. " words, more than APA's " ..
      kAbstractWords .. ". If the front matter of the first page no longer fits " ..
      "on it, it is scaled down to fit.")
  end
end

--- The whole of the body for a journal-mode pdf: the front matter, with the
--- author note at the foot of the first column, then the document.
function M.latex(meta, notes, corresponding, tail, blocks)
  warn_long_abstract(meta)
  local out = List:new {}
  out:insert(pandoc.Div(masthead_blocks(meta), pandoc.Attr("", { "JournalMasthead" })))
  local foot = footnote_blocks(meta, notes, corresponding)
  if #foot > 0 then
    out:insert(pandoc.Div(foot, pandoc.Attr("", { "JournalNote" })))
  end
  out:insert(raw([[\quasifirstpar]]))
  out:extend(tail)
  out:extend(blocks)
  return out
end

return M
