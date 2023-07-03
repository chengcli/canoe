The [canoe](https://github.com/chengcli/canoe) repository is public by default
but Github does not allow the creation of private forks for public repositories.

Thus, we recommend a public fork of [canoe](https://github.com/chengcli/canoe)
to create your own features. However, we understand that sometime a private fork is necessary.

The correct way of creating a private frok from a public repo is by duplicating the repo,
which is formally documented [here](https://help.github.com/articles/duplicating-a-repository/) by Github.

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

 1. Push the repository's Git Large File Storage objects to your mirror.
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

 1. Push and pull changes.
    When you push, do so on `origin` with `git push origin`.
    > Replace `<branch_naem>` with your the name of your development branch.
    ```bash
    git push origin <branch_name>
    ```
    When you want to pull changes from `upstream`, do so on `upstream` and `rebase` on top of your work.
    ```bash
    git fetch upstream
    git rebase upstream/master
    ```

 1. Finally, resolve the conflicts if any.
