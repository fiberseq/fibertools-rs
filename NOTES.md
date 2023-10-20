# How to install pytorch on conda with gpu

```
mamba install pytorch=1.12 torchvision torchaudio pytorch-cuda=11.6 -c pytorch -c nvidia
```

# How to make the docker image using the Dockerfile

```
docker build --progress=plain -f Dockerfile .
```
