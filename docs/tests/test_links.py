from pathlib import Path
import re
import itertools
from urllib import request

import pytest


@pytest.fixture
def docs_path():
    return Path('/home/alexandre/projects/neuro/neuro_pypes/docs/')


@pytest.fixture
def docs_contents(docs_path):
    doc_files = docs_path.glob('**/*.md')
    yield (doc.read_text() for doc in doc_files)


@pytest.fixture
def docs_urls(docs_contents):
    name_pattern = "[^]]+"
    url_pattern = "http[s]?://[^)]+"
    markup_pattern = '\[({0})]\(\s*({1})\s*\)'.format(name_pattern, url_pattern)
    http_regex = re.compile(markup_pattern, re.IGNORECASE)

    doc_urls = (http_regex.findall(content) for content in docs_contents)
    yield itertools.chain(*(url_group for url_group in doc_urls if url_group))


def test_urls(docs_urls):
    for url in docs_urls:
        _test_url(url[1])


def _test_url(url):
    print(url)
    response = request.urlopen(url)
    assert response.status == 200
