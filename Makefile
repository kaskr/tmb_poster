pdf: linreg rw
	pdflatex poster.tex

linreg:
	R CMD Sweave linreg.Rnw

rw:
	R CMD Sweave rw.Rnw

clean:
	rm -f *.so *.o *.pdf *~ *.out *.aux *.log rw.tex linreg.tex

