.DEFAULT_GOAL := analysis

requirements:
	pip install -r requirements.txt

makeref:
	python src/guides_to_ref.py

# process project's data
# with looper/pypiper/pipelines:
# see https://github.com/epigen/looper
# see https://github.com/epigen/pypiper
# see https://github.com/epigen/pipelines
process:
	looper run metadata/config.yaml

assign:
	python assign_gRNA_cells.py

collect:
	python src/collect_expression.py

analysis: assign collect
	python src/analysis.py

all: requirements makeref process assign collect analysis

.PHONY: requirements makeref process assign collect analysis all
