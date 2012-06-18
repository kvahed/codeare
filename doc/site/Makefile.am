# Generate local documentation?
if BUILD_DOCS

RESOURCES = codeare.css screen.css
PAGES = $(patsubst xml/%.xml,%.html,$(wildcard xml/*.xml))
%.html: xml/%.xml $(STYLESHEET)
	$(XSLTPROC) $(STYLESHEET) $< > $@

endif

EXTRA_DIST=$(PAGES) $(RESOURCES)

sitedir = $(datadir)/codeare/doc
site_DATA = $(EXTRA_DIST)

XSLTPROC = xsltproc --xinclude --nonet
STYLESHEET = xsl/codeare.org.xsl

# Generate local documentation?
if BUILD_DOCS

%.html: xml/%.xml
	$(XSLTPROC) $(STYLESHEET) $< > $@

endif
