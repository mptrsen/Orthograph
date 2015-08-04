doc: manual

manual:
	man -l -Tps doc/orthograph.man > doc/orthograph.ps
	ps2pdf doc/orthograph.ps doc/orthograph.pdf
	$(RM) doc/orthograph.ps
