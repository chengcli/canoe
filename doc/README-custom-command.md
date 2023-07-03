Git helps us to write our own custom commands. With these commands, you can write your own git tasks. These git commands can be written in any language you prefer like Ruby, JavaScript.

In this tutorial, I will explain how to write your own custom commands and I will give some examples of it. We will be using bash for this tutorial.

Writing custom git commands are simple as these three steps:

    Create a file with a name in the given format git-mycommand. And this file should be given executable access.

    The created file should be included in $PATH. This can be done in bashrc or in zshrc

    Then you can run the command git mycommand in any of the git repo you have.

Writing our first custom command

Let’s write a custom git command to add, commit, and then push the changes to remote with one command. This is just an example to write a custom command so don’t worry about the use-case of the command we are gonna build.
1. Creating a file for our command
It’s good to have a directory to contain all our custom commands. So that it will be organized. Let’s create a directory.

It’s good to have a directory to contain all our custom commands. So that it will be organized. Let’s create a directory.

mkdir my-git-custom-commands

Create a file named git-lazypush in the above created directory. And add the following code. Refer to the code in the gist here.

```bash
#!/bin/sh

message=$1 # First parameter will be the commit message
currentBranch=$(git symbolic-ref --short -q HEAD) # Getting the current branch

if [ ! -z "$1" ] # checking if the commit message is present. If not then aborting.
then
  git add .
  git commit -m "$message"
  git push origin $currentBranch
else
  echo "Commit message is not provided"
fi
```

Last but not least, make that file as executable

chmod +x git-lazypush

2. Add custom commands to $PATH

This is quite important. This will help Git to recognize the custom commands.

We are adding the directory to $PATH. Add the following line to bashrc or to zshrc

export PATH=$PATH:/your-directory/my-git-custom-commands

You can now source the file or start a new instance of the terminal

```bash
source ~/.zshrc
```

3. Running our custom command

With everything in place, it’s time to run our custom command. Go to any of your local git repo and run the following command

git lazypush “Your commit message“
