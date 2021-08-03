### main manuscript

all: gen_hmm.pdf DOapp/do_analysis.html CCapp/cc_analysis.html

R_OPTS=--no-save --no-restore --no-init-file --no-site-file # --vanilla, but without --no-environ

FIGS=Figs/fig1_genome_reconstr.pdf \
	 Figs/fig2_do_trmatrix.pdf \
	 Figs/fig3_do_qtl.pdf \
	 Figs/fig4_cc_xchr.pdf

DOapp/do_analysis.html: DOapp/do_analysis.Rmd
	cd $(<D);R $(R_OPTS) -e "rmarkdown::render('$(<F)')"

CCapp/cc_analysis.html: CCapp/cc_analysis.Rmd
	cd $(<D);R $(R_OPTS) -e "rmarkdown::render('$(<F)')"

Figs/%.pdf: R/%.R DOapp/do_analysis.html CCapp/cc_analysis.html
	cd $(<D);R CMD BATCH $(R_OPTS) $(<F)

# .tex to .pdf
gen_hmm.pdf: LaTeX/gen_hmm.tex gen_hmm.bib genetics.bst $(FIGS)
	cd LaTeX;pdflatex gen_hmm
	cd LaTeX;bibtex gen_hmm
	cd LaTeX;pdflatex gen_hmm
	cd LaTeX;pdflatex gen_hmm
	cd LaTeX;pdflatex gen_hmm
	\mv LaTeX/gen_hmm.pdf .

# Sweave to .tex
LaTeX/gen_hmm.tex: gen_hmm.Rnw DOapp/do_analysis.html CCapp/cc_analysis.html
	[ -d LaTeX ] || mkdir LaTeX
	[ -e LaTeX/genetics.bst ] || (cd LaTeX;ln -s ../genetics.bst)
	[ -e LaTeX/gen_hmm.bib ] || (cd LaTeX;ln -s ../gen_hmm.bib)
	[ -e LaTeX/Figs ] || (cd LaTeX;ln -s ../Figs)
	Rscript -e 'library(knitr); knit("gen_hmm.Rnw", "LaTeX/gen_hmm.tex")'
