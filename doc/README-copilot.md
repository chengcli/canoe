# Install nvim
```
yum install -y neovim python3-neovim
```

# Install node (21.x)
Checkout this website
[node distribution](https://github.com/nodesource/distributions#redhat-versions)
```
yum install https://rpm.nodesource.com/pub_21.x/nodistro/repo/nodesource-release-nodistro-1.noarch.rpm -y
sudo yum install nodejs -y --setopt=nodesource-nodejs.module_hotfixes=1
```

# Install Github Copilot
```
git clone https://github.com/github/copilot.vim.git \
  ~/.config/nvim/pack/github/start/copilot.vim
```
