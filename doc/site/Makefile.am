# Generate local documentation?
if BUILD_DOCS

RESOURCES = codeare.css screen.css sh_main.js sh_cpp.js sh_xml.js sh_style.css
PAGES = $(patsubst xml/%.xml,%.html,$(wildcard xml/*.xml))
%.html: xml/%.xml $(STYLESHEET)
	$(XSLTPROC) $(STYLESHEET) $< > $@

endif

EXTRA_DIST=$(PAGES) $(RESOURCES)

sitedir = $(DESTDIR)$(docdir)
site_DATA = $(EXTRA_DIST)

XSLTPROC = xsltproc --xinclude --nonet
STYLESHEET = xsl/codeare.org.xsl

# Generate local documentation?
if BUILD_DOCS

%.html: xml/%.xml
	$(XSLTPROC) $(STYLESHEET) $< > $@

endif
