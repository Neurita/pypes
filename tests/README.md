# Integration tests for neuro_pypes

## Installation

If you don't have the necessary software in your machine, you can clone neuro_ansible and launch the container:

```bash
git clone https://github.com/Neurita/neuro_ansible.git
git clone https://github.com/Neurita/multimodal_test_data

cd neuro_ansible

pipenv run packer build packer/neuro_docker-run.json

docker run -it --name neuro -v $PWD/../multimodal_test_data:/data neurita/neuro_docker:0.1 /bin/bash
```

## How to run the tests

```
pyenv activate neuro

```
