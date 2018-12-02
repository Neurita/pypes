
Build the neurita/neuro_docker container:

```
git clone https://github.com/Neurita/neuro_ansible.git
```

Install the dependencies of neuro_ansible. Then:

```
cd neuro_ansible
make docker-run
```
Exit from the container and delete it.

Run the container again with your options:
```
export DATA_DIR=$HOME/projects/multimodal_test_data
export PYPES_DIR=$HOME/projects/neuro_pypes

docker run -it -p 8888:8888 --name neuro -v $DATA_DIR:/data -v  $PYPES_DIR:/root/projects/neuro_pypes neurita/neuro_docker:0.2 /bin/bash
```

Inside the container:
```
pyenv activate neuro
pip install jupyter jupyterlab
```

```
jupyter lab --ip 0.0.0.0 --no-browser --allow-root
```
