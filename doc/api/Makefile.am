apidir = $(datadir)/codeare/doc/api
api_DATA = $(EXTRA_DIST)

static_docs = res/doxygen.css res/tabs.css
DOXY_FILE = doxy-codeare.conf
DOXY_STAMP = doxygen-build.stamp
doc_DATA = $(static_docs) $(DOXY_STAMP)
EXTRA_DIST = $(static_docs) $(DOXY_FILE).in

# Generate local API documentation?
if BUILD_API

doxygen-build.stamp: $(DOXY_FILE) $(static_docs) $(top_srcdir)/src/*.hpp $(top_srcdir)/src/matrix/*.hpp $(top_srcdir)/src/modules/*.hpp
	doxygen $(DOXY_FILE)
	touch $(DOXY_STAMP)
endif

install-data-local: doxygen-build.stamp
	$(mkinstalldirs) $(apidir)
	if test -d api; then \
	  for dir in `find api/* -type d`; do \
	    if test -d $$dir; then \
	      b=`basename $$dir`; \
	      $(mkinstalldirs) $(apidir)/$$b; \
	      for f in $$dir/*; do \
	        $(INSTALL_DATA) $$f $(apidir)/$$b; \
	      done \
	    fi \
	  done \
	fi

uninstall-local:
	rm -rf $(apidir)

clean-local:
	if test -d api; then rm -rf api; fi
	if test -f doxygen-build.stamp; then rm -f doxygen-build.stamp; fi
