## Writing custom git command

This article introduces two extremely useful custom commands, `git send` and `git done`.

Writing custom git commands are simple as these three steps:
    - Create a file with a name in the given format git-mycommand. And this file should be given executable access.
    - The created file should be included in $PATH. This can be done in bashrc or in zshrc
    - Then you can run the command git mycommand in any of the git repo you have.

Let’s write a custom git command to add, commit, and then push the changes to remote with one command.
Let's call this command `git send`, which writes a useless commit message "wip".

#### Step 1. Creating a file for our command
It’s good to have a directory to contain all our custom commands. So that it will be organized. Let’s create a directory.

```
mkdir .gitcustom
```

Create a file named `git-send` in the above created directory. And add the following code. Refer to the code in the gist here.

```bash
#!/bin/sh

currentBranch=$(git symbolic-ref --short -q HEAD) # Getting the current branch

git add .
git commit -m "wip"
git push origin $currentBranch
```

Last but not least, make that file as executable

```bash
chmod +x git-send
```

#### Step 2. Add custom commands to $PATH

This is quite important. This will help Git to recognize the custom commands.

We are adding the directory to $PATH. Add the following line to bashrc or to zshrc

```bash
export PATH=$PATH:<absolute_path>/.gitcustom
```
where `absolute_path` is the absolute path of the directory containing folder `.gitcustom`

You can now source the file or start a new instance of the terminal

```bash
source ~/.zshrc
```

#### Step 3. Running our custom command

With everything in place, it’s time to run our custom command. Go to any of your local git repo and run the following command

```bash
git send
```

### The git-done command

After your work has been merged into the main branch, you can safely delete
your working branch and rebase your local main branch. These steps can be
accomplished by the custom `git-done` command:

```bash
#!/bin/sh

currentBranch=$(git symbolic-ref --short -q HEAD) # Getting the current branch

git checkout main
git fetch
git rebase origin/main
git branch -D $currentBranch
```

Both commands are provided in the ``.gitcustom`` folder in ``Canoe``.
