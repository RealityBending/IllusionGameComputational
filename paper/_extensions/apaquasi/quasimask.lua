-- apaquasi: a masked document is set as apaquarto's plain manuscript.
--
-- The journal-style first page is made for a posted preprint: the doi, the
-- server, the date it went up, the links to its materials and the badges'
-- links to the osf all lead back to the authors, and a reviewer can look up a
-- page styled as a preprint. A double-blind submission is a manuscript, which
-- is what apaquarto's manuscript mode sets, masking included: the title alone
-- on the first page, then the abstract, then the text.
--
-- So mask: true turns journal mode into manuscript mode. Runs before
-- documentmode.lua, which is why it knows the long name for the mode too.

local stringify = pandoc.utils.stringify

function Meta(m)
  if not (m.mask and stringify(m.mask) == "true") then return nil end

  if m.documentmode then
    local mode = stringify(m.documentmode)
    if mode == "jou" or mode == "journal" then
      m.documentmode = "man"
    end
  end

  -- apaquarto prints the materials under the abstract whether masked or not,
  -- and a repository under a person's name gives them away. A view-only link
  -- made for the review is fine, so it is left in, and the writer told.
  if m["supplemental-materials"] and not m["suppress-supplemental-materials"] then
    quarto.log.warning("The document is masked, and its supplemental-materials " ..
      "are still printed under the abstract. Check that their links do not " ..
      "identify the authors (an anonymised view-only link, for example), or set " ..
      "suppress-supplemental-materials: true.")
  end
  return m
end
