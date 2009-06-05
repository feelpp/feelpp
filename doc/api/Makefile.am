# -*- makefile -*-
#
SUFFIXES	= .cpp .hpp .idl .c .h .f .F .o .moc


DOCSOURCES = $(shell find @top_srcdir@/../life  -name "*.[ch]pp") \
 $(shell find @srcdir@ -name "*.do[cx]")

EXTRA_DIST = Doxyfile \
 index.doc publis.doc  \
 $(foreach file,$(FIGS),fig/$(file)) $(PNGS)


CLEANFILES = .*.toc
DISTCLEANFILES = index-API.doc life.tag
MAINTAINERCLEANFILES = html

if WITH_DOXYGEN
%.png: fig/%.fig
	$(FIG2DEV) -Lpng $^ > $@

html: Doxyfile  $(PNGS) $(DOCSOURCES)
	: # Setting up
	touch index-API.doc .index.doc.toc
	test -e html || mkdir html
#	cp $(foreach file,$(PNGS),@srcdir@/$(file)) html/
	cp @srcdir@/life.css html/
	: # First pass
	$(DOXYGEN) Doxyfile >/dev/null
#	perl @top_srcdir@/tools/$(DOXYGEN)/index_create.pl life.tag index-API.doc
#	perl @top_srcdir@/tools/$(DOXYGEN)/toc_create.pl @srcdir@/index.doc @srcdir@/publis.doc
#	: # Second pass
#	$(DOXYGEN) Doxyfile >/dev/null
#	: # Post-processing
#	rm -f html/dir*
#	perl @top_srcdir@/tools/$(DOXYGEN)/$(DOXYGEN)_postprocessor.pl

#publish:
#	$(MAKE) -C .. sync-gforge
endif # WITH_DOXYGEN


install-data-local:
	@for file in `find html/ -type d` ; do \
	   echo "test -z \"$(DESTDIR)$(docdir)/$$file\" || $(mkdir_p) \"$(DESTDIR)$(docdir)/$$file\"";\
	   test -z "$(DESTDIR)$(docdir)/$$file" || $(mkdir_p) "$(DESTDIR)$(docdir)/$$file";\
	done
	@for file in `find html/ -type f` ; do \
	   echo " $(INSTALL_DATA) '$$file' '$(DESTDIR)$(docdir)/$$file'";\
	   $(INSTALL_DATA) $$file $(DESTDIR)$(docdir)/$$file ;\
	done

uninstall-local:
	for file in `find html/ -type f` ; do \
	   echo " rm -f '$(DESTDIR)$(docdir)/$$file'";\
	   rm -f $(DESTDIR)$(docdir)/$$file ;\
	done

.PHONY: html

include $(top_srcdir)/admin/dist-files.mk
