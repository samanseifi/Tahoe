# Build the technical report:  make   (requires latexmk + pdflatex)
LATEXMK := latexmk
LATEXMK_FLAGS := -pdf -interaction=nonstopmode -halt-on-error

.PHONY: all clean
all: technical_report.pdf

technical_report.pdf: technical_report.tex necking_validation/fig15_compare.png fig18_validation/fig18_x4_overlay.png
	$(LATEXMK) $(LATEXMK_FLAGS) technical_report.tex

clean:
	$(LATEXMK) -c technical_report.tex
