# Deploy an environment image

First, log out any existing Docker sessions:
```bash
<ctrl>+<D>
```

Use the following commands to check running containers:
```bash
make ps
```

Commit your changes to a new Docker image:
```bash
docker commit <container_id> ubuntu22.04-cuda12.9-py3.10-canoe:latest
```

Tag your image for a registry:
```bash
docker tag ubuntu22.04-cuda12.9-py3.10-canoe:latest docker.io/<DOCKERHUB_USER>/ubuntu22.04-cuda12.9-py3.10-canoe:YYYY-MM-DD
```

Log in and push your image to Docker Hub:
```bash
docker login
docker push docker.io/<DOCKERHUB_USER>/ubuntu22.04-cuda12.9-py3.10-canoe:YYYY-MM-DD
```

Use the image in manifest:
```
image: docker.io/<DOCKERHUB_USER>/ubuntu22.04-cuda12.9-py3.10-canoe:YYYY-MM-DD
imagePullPolicy: IfNotPresent
```
