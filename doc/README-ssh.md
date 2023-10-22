# configure login without password
1. Generate ssh keys
```
ssh-keygen
```

2. Copy ssh pub key to server
```
ssh-copy-id -i ~/.ssh/id_rsa.pub chengcli@ziggy.engin.umich.edu
```
