## Large file storage

This repo uses Git Large File (lfs) Storage to store benchmark data.
To install lfs for this repo, use
```
git lfs install
```
Sometimes, you will run into a authentication issue with git lfs.
Run the following command to check your environment:
```
git lfs env
```
If you find that the output looks like this
```
Endpoint=https://chengcli@github.com/chengcli/canoe.git/info/lfs (auth=none)
```
then you have encountered an authentication problem. Use the following step
to fix it:
```
git config --global lfs.https://github.com/chengcli/canoe.git/info/lfs.access basic
```
This command changes the autentication method from `none` to `basic`.

To pull files from lfs, use
```
git lfs pull
```

To check what files are stored in lfs, use
```
git lfs ls-files
```

To track large data file, use
```
git lfs track <file>
```


# Using git filter-repo

>Warning: If you run git filter-repo after stashing changes, you won't be able to retrieve your changes with other stash commands.
>Before running git filter-repo, we recommend unstashing any changes you've made.
>To unstash the last set of changes you've stashed, run `git stash show -p | git apply -R`.
>For more information, see [Git Tools - Stashing and Cleaning](https://git-scm.com/book/en/v2/Git-Tools-Stashing-and-Cleaning).

To illustrate how git filter-repo works, we'll show you how to remove your file with sensitive data from the history of your repository and add it to .gitignore to ensure that it is not accidentally re-committed.

  1. Install the latest release of the git filter-repo tool. You can install git-filter-repo manually or by using a package manager. For example, to install the tool with HomeBrew, use the brew install command.
  ```bash
  brew install git-filter-repo
  ```

  For more information, see INSTALL.md in the newren/git-filter-repo repository.

  2. If you don't already have a local copy of your repository with sensitive data in its history, clone the repository to your local computer.
  ```bash
  $ git clone https://github.com/YOUR-USERNAME/YOUR-REPOSITORY
  > Initialized empty Git repository in /Users/YOUR-FILE-PATH/YOUR-REPOSITORY/.git/
  > remote: Counting objects: 1301, done.
  > remote: Compressing objects: 100% (769/769), done.
  > remote: Total 1301 (delta 724), reused 910 (delta 522)
  > Receiving objects: 100% (1301/1301), 164.39 KiB, done.
  > Resolving deltas: 100% (724/724), done.
  ```

  3. Navigate into the repository's working directory.

  ```bash
  cd YOUR-REPOSITORY
  ```

  4. Run the following command, replacing PATH-TO-YOUR-FILE-WITH-SENSITIVE-DATA with the path to the file you want to remove, not just its filename. These arguments will:
  - Force Git to process, but not check out, the entire history of every branch and tag
  - Remove the specified file, as well as any empty commits generated as a result
  - Remove some configurations, such as the remote URL, stored in the .git/config file. You may want to back up this file in advance for restoration later.
  - Overwrite your existing tags

  ```bash
  $ git filter-repo --invert-paths --path PATH-TO-YOUR-FILE-WITH-SENSITIVE-DATA
  Parsed 197 commits
  New history written in 0.11 seconds; now repacking/cleaning...
  Repacking your repo and cleaning out old unneeded objects
  Enumerating objects: 210, done.
  Counting objects: 100% (210/210), done.
  Delta compression using up to 12 threads
  Compressing objects: 100% (127/127), done.
  Writing objects: 100% (210/210), done.
  Building bitmaps: 100% (48/48), done.
  Total 210 (delta 98), reused 144 (delta 75), pack-reused 0
  Completely finished after 0.64 seconds.
  ```
  > Note: If the file with sensitive data used to exist at any other paths (because it was moved or renamed), you must run this command on those paths, as well.

  5. Add your file with sensitive data to .gitignore to ensure that you don't accidentally commit it again.
  ```bash
  $ echo "YOUR-FILE-WITH-SENSITIVE-DATA" >> .gitignore
  $ git add .gitignore
  $ git commit -m "Add YOUR-FILE-WITH-SENSITIVE-DATA to .gitignore"
  > [main 051452f] Add YOUR-FILE-WITH-SENSITIVE-DATA to .gitignore
  >  1 files changed, 1 insertions(+), 0 deletions(-)
  ```

  Double-check that you've removed everything you wanted to from your repository's history, and that all of your branches are checked out.

  6. Once you're happy with the state of your repository, force-push your local changes to overwrite your repository on GitHub.com, as well as all the branches you've pushed up. A force push is required to remove sensitive data from your commit history.
  ```bash
  $ git push origin --force --all
  > Counting objects: 1074, done.
  > Delta compression using 2 threads.
  > Compressing objects: 100% (677/677), done.
  > Writing objects: 100% (1058/1058), 148.85 KiB, done.
  > Total 1058 (delta 590), reused 602 (delta 378)
  > To https://github.com/YOUR-USERNAME/YOUR-REPOSITORY.git
  >  + 48dc599...051452f main -> main (forced update)
  ```

  7. In order to remove the sensitive file from your tagged releases, you'll also need to force-push against your Git tags:
  ```bash
  $ git push origin --force --tags
  > Counting objects: 321, done.
  > Delta compression using up to 8 threads.
  > Compressing objects: 100% (166/166), done.
  > Writing objects: 100% (321/321), 331.74 KiB | 0 bytes/s, done.
  > Total 321 (delta 124), reused 269 (delta 108)
  > To https://github.com/YOUR-USERNAME/YOUR-REPOSITORY.git
  >  + 48dc599...051452f main -> main (forced update)
  ```

### Resources
[github](https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/removing-sensitive-data-from-a-repository)
