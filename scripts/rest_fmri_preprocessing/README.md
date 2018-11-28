
Run the neurita/neuro_docker container:

```
docker run -it -p 8888:8888 --name neuro -v $PWD/../multimodal_test_data:/data -v $PWD/../neuro_pypes:/root/projects/neuro_pypes neurita/neuro_docker:0.2 /bin/bash
```

Inside the container:
```
pyenv activate neuro
pip install jupyter jupyterlab
```

```
jupyter lab --ip 0.0.0.0 --no-browser --allow-root
```
