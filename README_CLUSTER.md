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
