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

#### Install kuberctl

1. Check this webpage for updates:
https://kubernetes.io/docs/tasks/tools/install-kubectl-linux/

2. Add kubernetes yum repository
```
cat <<EOF | sudo tee /etc/yum.repos.d/kubernetes.repo
[kubernetes]
name=Kubernetes
baseurl=https://pkgs.k8s.io/core:/stable:/v1.35/rpm/
enabled=1
gpgcheck=1
gpgkey=https://pkgs.k8s.io/core:/stable:/v1.35/rpm/repodata/repomd.xml.key
EOF
```

3. Install kubectl using yum
```
sudo yum install -y kubectl
```

#### Install k3s cluster (server)

1. Check this webpage for updates:
```
https://docs.k3s.io/quick-start
```

2. Download k3s and install
```
curl -sfL https://get.k3s.io | sh -
```

3. Copy kubeconfig to home directory
```
mkdir -p ~/.kube
sudo cp /etc/rancher/k3s/k3s.yaml ~/.kube/config
sudo chown $USER:$USER ~/.kube/config
chmod 600 ~/.kube/config
```

4. Check cluster info
```
kubectl cluster-info
```

5. Check node
```
kubectl get nodes
```

6. Find and copy node token
```
sudo cat /var/lib/rancher/k3s/server/node-token
```

7. Open network communiction ports
```
sudo firewall-cmd --permanent --add-port=6443/tcp
sudo firewall-cmd --permanent --add-port=8472/udp
sudo firewall-cmd --permanent --add-port=10250/tcp
sudo firewall-cmd --reload
sudo firewall-cmd --list-ports
```

8. (optional) uninstall k3s
```
k3s-killall.sh
k3s-uninstall.sh
```

#### Install k3s cluster (worker)
1. On any worker node, repeat the process of installing docker
2. Repeat the process of installing NVIDIA container tool kit

3. Verify server port
```
nc -vz dart9.engin.umich.edu 6443
```

4. Use the node token to install and join the server
```
curl -sfL https://get.k3s.io | K3S_URL=<SERVER_URL>:6443 K3S_TOKEN=<NODE_TOKEN> sh -
```
