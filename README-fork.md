The [canoe repository](https://github.com/chengcli/canoe) is public by default 
and Github does not allow the creation of private forks for public repositories.

We recommend a public fork of canoe to create your own features. However, sometime a
private fork is necessary.

The correct way of creating a private frok by duplicating the repo is documented [here](https://help.github.com/articles/duplicating-a-repository/).

For the canoe repo the commands are:

 1. Create a bare clone of the repository.
    (This is temporary and will be removed so just do it wherever.)
    ```bash
    git clone --bare https://github.com/chengcli/canoe.git
    ```

 2. [Create a new private repository on Github](https://help.github.com/articles/creating-a-new-repository/) 
    and name it `canoe-dev` (or anything else you would like to have).

 3. Mirror-push your bare clone to your new `canoe-dev` repository.
    > Replace `<username>` with your actual Github username in the url below.
    
    ```bash
    cd canoe.git
    git push --mirror https://<username>@github.com/<username>/canoe-dev.git
    ```

 4. Remove the temporary local repository you created in step 1.
    ```bash
    cd ..
    rm -rf canoe.git
    ```
    
 5. You can now clone your `canoe-dev` repository on your machine
    ```bash
    git clone https://<username>@github.com/<username>/canoe-dev.git
    ```
    Note that `git lfs` objects cannot not be cloned due to permission issues.
    But they are inconsequential.
   
 6. Add the original repo as remote to fetch future changes.
    ```bash
    git remote add upstream https://github.com/chengcli/canoe.git
    ```
    Now, you can list all your remotes with `git remote -v`. You should see:
    ```
    origin	https://<username>@github.com/<username>/canoe-dev (fetch)
    origin	https://<username>@github.com/<username>/canoe-dev (push)
    upstream	https://github.com/chengcli/canoe.git (fetch)
    upstream	https://github.com/chengcli/canoe.git (push)
    ```
    > When you push, do so on `origin` with `git push origin`.
   
    > When you want to pull changes from `upstream` you can just fetch the remote and rebase on top of your work.
    ```bash
      git fetch upstream
      git rebase upstream/master
      ```
      And solve the conflicts if any
