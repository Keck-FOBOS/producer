# Makefile for Sphinx documentation
#

# You can set these variables from the command line.
SPHINXAPI     = sphinx-apidoc
SPHINXAPIOPT  = members,private-members,undoc-members,show-inheritance
SPHINXOPTS    = -aE -w ./sphinx_warnings.out
SPHINXBUILD   = sphinx-build
BUILDDIR      = _build
STATICDIR     = _static
DOCTREE       = $(BUILDDIR)/doctrees
LOCALFILES    = $(BUILDDIR)/* api/*.rst sphinx_warnings.out

# Internal variables.
ALLSPHINXOPTS   = -d $(DOCTREE) $(SPHINXOPTS) .

clean:
	rm -rf $(LOCALFILES)

apirst:
	SPHINX_APIDOC_OPTIONS=$(SPHINXAPIOPT) $(SPHINXAPI) --separate -o ./api ../producer ../producer/tests/* ../producer/scripts/* 
	python ./scripts/write_script_help.py
	python ./scripts/fobos_layout_table.py

htmlonly:
	$(SPHINXBUILD) -b html $(ALLSPHINXOPTS) $(BUILDDIR)/html

html:
	make apirst
	make htmlonly


