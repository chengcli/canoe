# Set up a kubernetes multi-node GPU cluster

This tutorial walks you through setting up a **multi-node GPU cluster** using k3s,
enabling you to scale beyond a single machine and expand your computational capabilities.

It consolidates guidance from multiple sources. Some references may evolve over time, 
so if you encounter issues, be sure to consult the latest versions of the documentation below:

1. [Docker](https://docs.docker.com/engine/install/rhel/)
2. [NVIDIA Container Toolkit](https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/latest/install-guide.html)
3. [K3s cluster](https://docs.k3s.io/quick-start)
4. [NVIDIA Device Plugin](https://github.com/NVIDIA/k8s-device-plugin?tab=readme-ov-file)
5. [Tutorial to set up a single node cluster](https://www.radicalgeek.co.uk/adding-a-gpu-node-to-a-k3s-cluster/)

**Note:** The [kind](https://kind.sigs.k8s.io/docs/user/quick-start#installing-with-a-package-manager)
cluster supports **single-node deployments** only and is therefore not suitable for multi-node GPU setups.

## Install docker on Redhat

1. Install the `dnf-plugins-core` package:
```
sudo dnf -y install dnf-plugins-core
sudo dnf config-manager --add-repo https://download.docker.com/linux/rhel/docker-ce.repo
```

2. Install the lagest version
```
sudo dnf install docker-ce docker-ce-cli containerd.io docker-buildx-plugin docker-compose-plugin
```

3. Start docker engine
```
sudo systemctl enable --now docker
```

4. Verify that docker has been successfully installed
```
sudo docker run hello-world
```

5. Add user to the docker group
```
sudo usermod -aG docker $USER
```

6. Log out and log back in to take effect and validate with
```
docker ps
```

You should see the following without errors:
```
CONTAINER ID   IMAGE     COMMAND   CREATED   STATUS    PORTS     NAMES
```

## Install NVIDIA Container Toolkit

1. Configure the production repository
```
curl -s -L https://nvidia.github.io/libnvidia-container/stable/rpm/nvidia-container-toolkit.repo | \
  sudo tee /etc/yum.repos.d/nvidia-container-toolkit.repo
```

2. Install the NVIDIA Container Toolkit packages
```
export NVIDIA_CONTAINER_TOOLKIT_VERSION=1.18.2-1
  sudo dnf install -y \
      nvidia-container-toolkit-${NVIDIA_CONTAINER_TOOLKIT_VERSION} \
      nvidia-container-toolkit-base-${NVIDIA_CONTAINER_TOOLKIT_VERSION} \
      libnvidia-container-tools-${NVIDIA_CONTAINER_TOOLKIT_VERSION} \
      libnvidia-container1-${NVIDIA_CONTAINER_TOOLKIT_VERSION}
```

3. Let docker use NVIDIA's **containerd** runtime
```
sudo nvidia-ctk runtime configure --runtime=containerd
```

4. Restart containerd servier
```
sudo systemctl restart containerd
```

## Pull NVIDIA docker images

1. Pull docker images with NVIDIA CUDA
```
docker pull nvidia/cuda:12.8.0-devel-ubuntu22.04
```

2. Check GPUs are correctly recognized by docker
```
docker run --rm --gpus all nvidia/cuda:12.8.0-runtime-ubuntu22.04 nvidia-smi
```

## Install k3s cluster (server)

1. Download k3s and install
```
curl -sfL https://get.k3s.io | sh -
```

2. Copy kubeconfig to home directory
```
mkdir -p ~/.kube
sudo cp /etc/rancher/k3s/k3s.yaml ~/.kube/config
sudo chown $USER:$USER ~/.kube/config
chmod 600 ~/.kube/config
```

3. Check cluster info
```
kubectl cluster-info
```

4. Check node
```
kubectl get nodes
```

5. Find and copy node token
```
sudo cat /var/lib/rancher/k3s/server/node-token
```

6. Open network communiction ports
```
sudo firewall-cmd --permanent --add-port=6443/tcp
sudo firewall-cmd --permanent --add-port=8472/udp
sudo firewall-cmd --permanent --add-port=10250/tcp
sudo firewall-cmd --reload
sudo firewall-cmd --list-ports
```

7. (Optional) uninstall k3s
```
k3s-killall.sh
k3s-uninstall.sh
```

## Join k3s cluster (worker)

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

## Let k3s recognize your GPU resource

1. Check cluster setup
```
kubectl get nodes
```

2. Describe a node
```
kubectl describe node csrwks2024-0242.engin.umich.edu
```
GPU resource is not allocatable now.

3. Label roles
```
kubectl label node csrwks2024-0242.engin.umich.edu node-type=worker
kubectl label node csrwks2024-0243.engin.umich.edu node-type=worker
kubectl label node csrwks2024-0244.engin.umich.edu node-type=worker
```
4. Create an nvidia runtime class
```
cat > runtime-class.yaml <<'EOF'
apiVersion: node.k8s.io/v1
kind: RuntimeClass
metadata:
  name: nvidia
handler: nvidia
EOC
```

5. Deploy to the cluster
```
kubeclt apply -f nvidia-runtimeclass.yaml
```

6. Install NVIDIA-PLUGIN
```
kubectl create -f nvidia-device-plugin.yaml
```

7. Check plugin deamon running
```
kubectl get pods -n kube-system | grep nvidia
```

8. Check GPU resource
```
kubectl get nodes -o custom-columns=NAME:.metadata.name,GPU_CAPACITY:.status.capacity.'nvidia\.com/gpu',GPU_ALLOCATABLE:.status.allocatable.'nvidia\.com/gpu'
```

9. (Optional) Remove daemon set
```
kubectl delete daemonset nvidia-device-plugin-daemonset -n kube-system
```

## Final test

1. Create a test file:
```
cat > gpu-test.yaml <<'EOF'
apiVersion: v1
kind: Pod
metadata:
  name: gpu-test
spec:
  runtimeClassName: nvidia
  containers:
  - name: cuda-test
    image: nvidia/cuda:12.3.2-base-ubuntu22.04
    command: ["nvidia-smi"]
    resources:
      limits:
        nvidia.com/gpu: 1
  restartPolicy: OnFailure
EOF
```

2. Create the pod
```
kubectl apply -f gpu-test.yaml
```

3. Watch it run
```
kubectl get pod gpu-test -w
```

4. Get actual result
```
kubectl logs gpu-test
```

5. Success would look like
```
+-----------------------------------------------------------------------------+
| NVIDIA-SMI 570.xx.xx    Driver Version: 570.xx    CUDA Version: 12.x        |
|-------------------------------+----------------------+----------------------+
| GPU  Name        Persistence-M| Bus-Id        Disp.A | Volatile Uncorr. ECC |
| Fan  Temp  Perf  Pwr:Usage/Cap|         Memory-Usage | GPU-Util  Compute M. |
|                               |                      |               MIG M. |
|===============================+======================+======================|
|  0  RTX A6000 / A100 / etc...                                      |
+-----------------------------------------------------------------------------+
```
