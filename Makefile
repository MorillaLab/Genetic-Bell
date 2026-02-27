.PHONY: install lint model em_peaks clean help

help:
	@echo "Genetic-Bell — available commands:"
	@echo "  make install    Install R package dependencies"
	@echo "  make lint       Lint R source files"
	@echo "  make model      Run the core Genetic-Bell model"
	@echo "  make em_peaks   Run EM peak detection"
	@echo "  make clean      Remove cache and temporary files"

install:
	Rscript requirements.R

lint:
	Rscript -e "if (!requireNamespace('lintr', quietly=TRUE)) install.packages('lintr'); \
	            lintr::lint_dir('model'); lintr::lint_dir('EM_peaks')"

model:
	Rscript -e "source('model/genetic_bell_model.R')"

em_peaks:
	Rscript -e "source('EM_peaks/em_peak_detection.R')"

clean:
	find . -name "*.Rhistory" -delete 2>/dev/null; true
	find . -name ".RData" -delete 2>/dev/null; true
	find . -name "Rplots.pdf" -delete 2>/dev/null; true
	find . -name ".DS_Store" -delete 2>/dev/null; true
