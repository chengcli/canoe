# How to private fork [canoe](https://github.com/chengcli/canoe)
The [canoe](https://github.com/chengcli/canoe) repository is public by default
but Github does not allow the creation of private forks for public repositories.

Thus, we recommend a public fork of [canoe](https://github.com/chengcli/canoe)
to create your own features. However, we understand that sometime a private fork is necessary.

The correct way of creating a private frok from a public repo is by duplicating the repo,
which is formally documented [here](https://help.github.com/articles/duplicating-a-repository/) by Github.

Make sure that you have at least git 2.39 and git lfs 3.2.
Use `git --version` and `git lfs --version` to check your version. The method below does
not work for lower version of git and git lfs.

For the [canoe](https://github.com/chengcli/canoe) repo, the specific commands are:

 1. Create a bare clone of [canoe](https://github.com/chengcli/canoe).
    (This is temporary and will be removed so just do it wherever)
    ```bash
    git clone --bare https://github.com/chengcli/canoe.git
    ```

 1. [Create a new private repository on Github](https://help.github.com/articles/creating-a-new-repository/)
    and name it `canoe-dev`, or anything else you would like to have.

 1. Pull in the repository's Git Large File Storage objects.
    ```bash
    cd canoe.git
    git lfs fetch --all
    ```

 1. Mirror-push your bare clone to your new `canoe-dev` repository.
    > Replace `<username>` with your actual Github username in the url below.

    ```bash
    git push --mirror https://<username>@github.com/<username>/canoe-dev.git
    ```

 1. Push the repository's Git Large File Storage objects to your fork.
    ```bash
    git lfs push --all https://<username>@github.com/<username>/canoe-dev.git
    ```

 1. Remove the temporary local repository you created in step 1.
    ```bash
    cd ..
    rm -rf canoe.git
    ```

 1. You can now clone your `canoe-dev` repository to your machine
    ```bash
    git clone https://<username>@github.com/<username>/canoe-dev.git
    ```

 1. Add the original repo as `upstream` to fetch future changes.
    ```bash
    git remote add upstream https://github.com/chengcli/canoe.git
    ```

 1. Now, you can list all your remotes with `git remote -v`. You should see:
    ```bash
    origin	https://<username>@github.com/<username>/canoe-dev (fetch)
    origin	https://<username>@github.com/<username>/canoe-dev (push)
    upstream	https://github.com/chengcli/canoe.git (fetch)
    upstream	https://github.com/chengcli/canoe.git (push)
    ```

 1. When you push, do so on `origin` with `git push origin`:
    > Replace `<branch_naem>` with your the name of your development branch.
    ```bash
    git push origin <branch_name>
    ```

 1. Finally, resolve the conflicts if any.

# How to fetch updates from upstream

If `upstream` has updates without lfs object, you can simply pull changes
from `upstream` by:
```bash
git fetch upstream
git rebase upstream/master
```
When `upstream` has new lfs object, do these extra steps before fetch:

 1. Create another bare clone of [canoe](https://github.com/chengcli/canoe).
     ```bash
     git clone --bare https://github.com/chengcli/canoe.git
     ```

 1. Pull in the updated Git Large File Storage objects.
    ```bash
    cd canoe.git
    git lfs fetch --all
    ```

 1. Push the repository's Git Large File Storage objects to your fork.
    ```bash
    git lfs push --all https://<username>@github.com/<username>/canoe-dev.git
    ```
Then, `fetch` and `rebase` will bring in the updates from `upstream`
```bash
git fetch upstream
git rebase upstream/master
```

# How to resolve conflicts
If `rebase` fails because of conflicting changes, you have four options:
 1. accept the changes from `upstream`
 1. accept local changes
 1. resolve conflicts manually
 1. revert to the state before rebase

## Over changes from one-side
Here, you would need the concept of `ours` and `theirs`.

- To accept the changes from `upstream` (if you know that `upstream` has the most updated code), use
  ```bash
  git rebase -Xours upstream/main
  ```
This command will override the local changes with updates from `upstream`.

- To accept local changes, use
  ```bash
  git rebase -Xtheirs upstream/main
  ```
This command will override the changes from `upstream` and with updates from local.

Note that during `rebase`, `ours` and `theirs` may appear swapped;
`ours` gives the version from the branch the changes are rebased onto, while
`theirs` gives the version from the branch that holds your work that is being rebased.

This is because `rebase` is used in a workflow that treats the history at the remote as the shared canonical one,
and treats the work done on the branch you are rebasing as the third-party work to be integrated,
and you are temporarily assuming the role of the keeper of the canonical history during the rebase.

As the keeper of the canonical history, you need to view the history from the remote as `ours` (i.e. "our shared canonical history"),
while what you did on your side branch as `theirs` (i.e. "one contributor’s work on top of it").

## Abort rebasing
You may abort unfinished rebasing by:
```bash
git rebase --abort
```
This command rolls back your repository to the state before `rebase`.

## Manually resolve the conflicts
TO BE FILLED
