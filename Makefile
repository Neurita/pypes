.PHONY: help clean clean-pyc clean-build list test test-dbg test-cov test-all coverage release sdist install deps develop upload path minor major

project-name = neuro_pypes

version-var := "__version__ = "
version-string := $(shell grep $(version-var) $(project-name)/version.py)
version := $(subst __version__ = ,,$(version-string))

help:
	@echo "clean-build - remove build artifacts"
	@echo "clean-pyc - remove Python file artifacts"
	@echo "lint - check style with flake8"
	@echo "test - run tests quickly with the default Python"
	@echo "test-cov - run tests with the default Python and report coverage"
	@echo "test-dbg - run tests and debug with pdb"
	@echo "testall - run tests on every Python version with tox"
	@echo "coverage - check code coverage quickly with the default Python"
	@echo "docs - generate Sphinx HTML documentation"
	@echo "install - install"
	@echo "develop - install in development mode"
	@echo "deps - install dependencies"
	@echo "dev_deps - install dependencies for development"
	@echo "release - package a release in wheel and tarball and upload it to PyPI"
	@echo "patch - bumpversion patch"
	@echo "minor - bumpversion minor"
	@echo "major - bumpversion major"

install:
	pipenv run python setup.py install
	pipenv install

develop:
	pipenv run python setup.py develop
	pipenv install --dev --skip-lock

clean: clean-build clean-pyc clean-pyenv

clean-build:
	rm -fr build/
	rm -fr dist/
	rm -fr *.egg-info

clean-pyc:
	find . -name '*.pyc' -exec rm -f {} +
	find . -name '*.pyo' -exec rm -f {} +
	find . -name '*~' -exec rm -f {} +
	find . -name '__pycache__' -exec rm -rf {} +
	find . -name '*.log*' -delete

clean-pyenv:
	pipenv --rm

lint:
	pipenv run flake8 $(project-name)/

test:
	pipenv run py.test -v

test-cov:
	pipenv run py.test --cov-report term-missing --cov=$(project-name)

test-dbg:
	pipenv run py.test --pdb

test-all:
	pipenv run tox

coverage:
	pipenv run coverage run --source $(project-name) setup.py test
	pipenv run coverage report -m

docs:
	rm -f docs/$(project-name).rst
	rm -f docs/modules.rst
	sphinx-apidoc -o docs/ $(project-name)
	$(MAKE) -C docs clean
	$(MAKE) -C docs html
	open docs/_build/html/index.html

tag: clean
	@echo "Creating git tag v$(version)"
	git tag v$(version)
	git push --tags

log:
	pipenv run gitchangelog

patch:
	pipenv run bumpversion patch

minor:
	pipenv run bumpversion minor

major:
	pipenv run bumpversion major

build:
	pipenv run python setup.py sdist --formats gztar bdist_wheel

upload:
	pipenv run python setup.py sdist upload

release: clean-build clean-pyc build upload tag
