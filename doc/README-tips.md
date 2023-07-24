
## Quick tips
- undo a "git add"
```
git reset <file>
```
- undo change to a file
```
git restore <file>
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
