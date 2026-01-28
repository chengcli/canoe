## This tutorial helps you setting up a kubernetes GPU cluster

#### Install docker on Redhat

1. Checkout this webpage for update:
https://docs.docker.com/engine/install/rhel/

2. Install the `dnf-plugins-core` package:
```
sudo dnf -y install dnf-plugins-core
sudo dnf config-manager --add-repo https://download.docker.com/linux/rhel/docker-ce.repo
```

3. Install the lagest version
```
sudo dnf install docker-ce docker-ce-cli containerd.io docker-buildx-plugin docker-compose-plugin
```

4. Start docker engine
```
sudo systemctl enable --now docker
```

5. Verify that docker has been successfully installed
```
sudo docker run hello-world
```

6. Add user to the docker group
```
sudo usermod -aG docker $USER
```

7. Log out and log back in to take effect and validate with
```
docker ps
```

You should not see the following without errors:
```
CONTAINER ID   IMAGE     COMMAND   CREATED   STATUS    PORTS     NAMES
```

#### Install NVIDIA Container Toolkit

1. Checkout this webpage for update:
```
https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/latest/install-guide.html
```

2. Configure the production repository
```
curl -s -L https://nvidia.github.io/libnvidia-container/stable/rpm/nvidia-container-toolkit.repo | \
  sudo tee /etc/yum.repos.d/nvidia-container-toolkit.repo
```

3. Install the NVIDIA Container Toolkit packages
```
export NVIDIA_CONTAINER_TOOLKIT_VERSION=1.18.2-1
  sudo dnf install -y \
      nvidia-container-toolkit-${NVIDIA_CONTAINER_TOOLKIT_VERSION} \
      nvidia-container-toolkit-base-${NVIDIA_CONTAINER_TOOLKIT_VERSION} \
      libnvidia-container-tools-${NVIDIA_CONTAINER_TOOLKIT_VERSION} \
      libnvidia-container1-${NVIDIA_CONTAINER_TOOLKIT_VERSION}
```

4. Let docker use NVIDIA's container runtime
```
sudo nvidia-ctk runtime configure --runtime=docker
```

5. Restart docker servier
```
sudo systemctl restart docker
```

#### Pull NVIDIA docker images

1. Pull docker images:
```
docker pull nvidia/cuda:12.8.0-devel-ubuntu22.04
```

2. Check GPUs are correctly recognized by docker
```
docker run --rm --gpus all nvidia/cuda:12.8.0-runtime-ubuntu22.04 nvidia-smi
```
