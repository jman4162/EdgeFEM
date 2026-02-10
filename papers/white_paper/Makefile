# Makefile for EdgeFEM Technical White Paper
#
# Usage:
#   make          - Build the PDF
#   make clean    - Remove auxiliary files
#   make distclean - Remove all generated files including PDF

MAIN = edgefem_methodology
LATEX = pdflatex
BIBTEX = bibtex

# LaTeX flags for better error reporting
LATEXFLAGS = -interaction=nonstopmode -halt-on-error

.PHONY: all clean distclean figures

all: $(MAIN).pdf

$(MAIN).pdf: $(MAIN).tex bibliography.bib
	$(LATEX) $(LATEXFLAGS) $(MAIN)
	$(BIBTEX) $(MAIN) || true
	$(LATEX) $(LATEXFLAGS) $(MAIN)
	$(LATEX) $(LATEXFLAGS) $(MAIN)
	@echo "Build complete: $(MAIN).pdf"

# Quick build without bibliography (for drafting)
quick:
	$(LATEX) $(LATEXFLAGS) $(MAIN)
	@echo "Quick build complete: $(MAIN).pdf"

# Generate figures from test data
figures:
	MPLBACKEND=Agg python3 generate_figures.py
	@echo "Figures generated in figures/"

clean:
	rm -f $(MAIN).aux $(MAIN).bbl $(MAIN).blg $(MAIN).log $(MAIN).out
	rm -f $(MAIN).toc $(MAIN).lof $(MAIN).lot $(MAIN).fls $(MAIN).fdb_latexmk
	rm -f $(MAIN).synctex.gz
	@echo "Auxiliary files cleaned"

distclean: clean
	rm -f $(MAIN).pdf
	@echo "All generated files cleaned"

# Watch for changes and rebuild (requires fswatch or inotifywait)
watch:
	@echo "Watching for changes to $(MAIN).tex..."
	@if command -v fswatch >/dev/null 2>&1; then \
		fswatch -o $(MAIN).tex bibliography.bib | xargs -n1 -I{} make; \
	elif command -v inotifywait >/dev/null 2>&1; then \
		while true; do \
			inotifywait -q -e modify $(MAIN).tex bibliography.bib; \
			make; \
		done; \
	else \
		echo "Install fswatch (macOS) or inotify-tools (Linux) for watch mode"; \
	fi

# Word count (approximate)
count:
	@echo "Approximate word count:"
	@detex $(MAIN).tex 2>/dev/null | wc -w || \
		echo "Install detex for word count, or use: texcount $(MAIN).tex"

# Check for common LaTeX issues
check:
	@echo "Checking for common issues..."
	@grep -n "\\\\cite{" $(MAIN).tex | head -20 || true
	@echo ""
	@echo "Undefined references (if any):"
	@grep -i "undefined" $(MAIN).log 2>/dev/null || echo "None found"
	@echo ""
	@echo "Missing citations (if any):"
	@grep -i "citation" $(MAIN).log 2>/dev/null | grep -i "undefined" || echo "None found"
