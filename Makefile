.PHONY: all dist doc test clean install uninstall

all: doc dist

docs/pdf/usage.pdf: docs/usage.rst
	pandoc --toc docs/usage.rst -o docs/pdf/usage.pdf


docs/pdf/method.pdf: docs/method.rst
	pandoc --toc docs/method.rst -o docs/pdf/method.pdf 

doc: docs/pdf/usage.pdf docs/pdf/method.pdf

clean-doc:
	rm docs/pdf/*

dist:
	python setup.py sdist

clean-dist:
	rm dist/*

clean: clean-doc clean-dist

install:
	pip install dist/tefingerprint-*.tar.gz

uninstall:
	pip uninstall -y tefingerprint

test: 
	python -m pytest -v ./

