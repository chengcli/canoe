
## Quick tips
- undo a "git add"
```
git reset <file>
```
- undo change to a file
```
git restore <file>
```
- undo a "git rebase"
```
git reset --hard ORIG_HEAD
```
- hard reset to remote
```
git reset --hard origin/main
```
- aliasing a git command
```
git config --global alias.co checkout
```
- clean *ALL* untracked files (this will cause unrecoverable damage, make sure what you
  are doing.
```
git clean -fd
```
- track a remote branch
Older versions of git:
```
git checkout --track <origin/branch>
```
Newer versions of git:
```
git checkout <branch>
```
- find where your package is installed
```
brew --prefix <package>
```
- cache your github password
```
git config --global credential.helper cache
```
- change the default password cache timeout,
```
git config --global credential.helper 'cache --timeout=86400'
```
- delete a remote branch
```
git push origin --delete <branch>
```
- configure your git signature
```
git config --global user.name "stormy/cli"
git config --global user.email "chengcli@umich.edu"
```
