# Canoe: Comprehensive Atmosphere N' Ocean Engine

[![build](https://github.com/chengcli/canoe/actions/workflows/main.yml/badge.svg)](https://github.com/chengcli/canoe/actions/workflows/main.yml)
[![build](https://github.com/chengcli/canoe/actions/workflows/mac.yaml/badge.svg)](https://github.com/chengcli/canoe/actions/workflows/mac.yaml)
[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

## Install system libraries and toolchain
Canoe can be installed on either a Linux distribution or on MacOS. Open a Linux or Mac terminal,
you can clone this repo using the following command:
```
git clone https://github.com/chengcli/canoe
```
This will copy all source files into your local computer. You will need to install a few
system libraries before installing canoe. All following instructions are executed under
the `canoe/` directory, which is referred to as the `root`.

### MacOS Installation Guide
We assume that [homebrew](https://brew.sh/) is already installed on your Mac and we will
use `brew` to install required system libraries. The system libraries are listed in
`Brewfile` at the `root`. To install them all, execute
```
brew bundle
```
### Ubuntu Linux Installation Guide
On a Ubuntu linux system, use `apt` to install
```
sudo apt install clang-format cmake
```

### Redhat Linux Installation Guide
On a Redhat linux system, use `yum` to install
```
sudo yum install clang-tools-extra cmake
```

## Install python libraries
All needed python libraries are collected in `requirements.txt`. We suggest using a
python [virtual environment](https://docs.python.org/3/library/venv.html) to install
these packages. If you are already using a virtual enviroment, install python packages
by
```
pip install -r requirements.txt
```
Otherwise, to create a python virtual environment:
```
python -m venv pyenv
```
This command will create an environment named `pyenv` in your current directory. Then, you
can use the previous command to install the python packages.

## Install pre-commit
Register your `pre-commit` hooks using
```
pre-commit install
```
The contributor's guide explains the meaning of `pre-commit`.

## How to build and test
After you completed the installation steps, you can build the canoe library.
The easiest way is to build it in-place, meaning that the build (binary files) are
located under `root`. To do so, make a new directory named `build`
```
mkdir build
```
All build files will be generated and placed under this directory. It is completely safe
to delete the whole directory if you want another build. `cd` to build and `cmake`

```
cd build
cmake ..
```
This command tells the cmake command to look for `CMakeFiles.txt` in the parent directory,
and start configuring the compile environment. Then compile the code by
```
make -j4
```
This comman will use 4 cores to compile the code in parallel. Once complete, all executable
files will be placed in `build/bin`.

## Contributors' Guide
We recommend the following naming conventions and workflow for all developers of this repo.
This ensures that the contribution can be successfully integrated into the existing code base
without breaking others' functionalities. We have configured this repo such that your contribution
will pass through several automated tests using [pre-commit](https://pre-commit.com/) hooks and
language lints such as [cpplint](https://github.com/cpplint/cpplint).

`pre-commit` hooks are triggered when you perform `git commit`.
The configuration file is located at [.pre-commit-config.yaml](https://github.com/chengcli/canoe/blob/main/.pre-commit-config.yaml).
The hooks are set of rules to format code automatically so that your code looks nice and
tidy. The hooks will perform the changes for your and you will need to add the files
changed the hooks again to `git` using `git add .`. Normally, you do not need to check
the files changed by the hooks.

Additionally, language lints checks your language style and gives suggestions. We
suggest fixing your code according to the suggestions exactly. Your code won't be able
to merge into the main branch without passing the lints. Here are a few naming
conventions to let you pass lints easier and help you write better and readable codes.

### Name your folders and files
- use a singular simple noun for folder names
- avoid compound nouns for folder names
- you can use compound nouns or phrases for files names
- individual words in a compound nouns for phrase should be concatenated by underscore

### Name your variables and classes
- use low case letters for variables
- you can use compound nouns or phrases
- compound nouns for phrase should be concatenated by underscore. This is called the *snake case*
- if a compound noun is used for class name, capitalize each word in the compound noun.
  This is called the *upper camel case*
- if a variable is a private member of a class, append the variable name with an
  underscore

### Name your funcions
- functions names are usualy verbal phrases
- standalone functions should use *snake case*
- public class member functions should use *upper camel case*
- private class member functions should use *lower camel case* in which the first word
  in a phase is not capitalized and the rest words are capitalized

### From C to C++ coding style
- use full path name starting from the `src` folder for include guard
- use `snprintf` instead of `sprintf`
- use `rand_r` instead of `rand`
- use `strtok_r` instead of `strtok`

## Git workflow
We adopt the idea of **linear history** and **squash merging**, meaning that there will
only be one permanent branch (main) and only way to push to main is by submitting a Pull
Request (PR). The main branch is protected such that no one can push to main directly.
**Linear history** ensures that the main branch is clean and tidy. **Sqush merging**
means that **the smallest unit of change is a PR, not a commit**.

This may be different from many people's individual workflow where the smallest unit is
a commit. For a collaborative project, commits are too fine-grainded and do not track
the issues very well. We want to make sure that each stage in the history solves a problem
that people can trace back and understand the context of this problem. Another way to
say it, the development is *issue driven*. The git workflow goes as the following:

### 1. Submit an issue ticket on Github website
Before you perform any work, think about what issue you want to solve by yourself or you
want other developers to solve it for you. When you have an idea of what this issue is
about, you create an issue ticket on Github. If you are going to solve the issue by
yourself, you can simply write a brief tile and not worrying about the details. If you
want some other developer to solve it for you, you may want to elaborate this issue for
others.

### 2. Create a new branch locally
If you have created or have seen an issue and you want to solve it by writing some new
code, you begin with creating a new branch locally.
```
git checkout -b <branch_name>
```
`branch_name` should begin with your `username` followed by a slash `/` followed by the
nature of this issue `job`. For example, `cli/add_cpptest_case` is a good branch name
indicating that this branch is created by user `cli` and works on adding cpp test cases.

### 3. Work on this branch
You can work on this branch by creating a document file under `doc/` folder.
This document should begin with a paragraph stating what is the problem and how you are
going to solve it. As you make progress, keep updating the document.

### 4. Update the `.gitignore` file
The `.gitignore` file helps you keep your working directory clean.
Each folder can have its own `.gitignore` file. This file records the files that
you don't wish to be added to the git system. For example, model output files should not
be added to the git. Ideally, when you execute `git status`, there should be no untracked
files in your git working directory.

### 5. Add changed files to git system
When you want to take a pause on working on the issue, you should add your changes to
the git system. This is done by using:
```
git add .
```
This command adds all untracked files to git excluding the files listed in the `.gitignore` file.

### 6. Commit your message
After you add the changed files, you can use
```
git status
```
to take a look at the changed files. If you accidentally add filles that you do not wish to add.
There is a way to *undo* add. Look for **Quick tips** section in this README file.
Then use the following command to commit your changes locally

```
git commit -m <message>
```
The content of `message` is not important. All commits within a PR (pull request) will
be squashed later into one message, you write later. You can either use a meaning
message like "work on XXX", "working XXX", or an unmeanding message like "wip", which
stands for "work in progress".

### 7. Upload your branch to GitHub
The previous command only commits the changes to your local machine. You can push your
changes to Github using
```
git push origin <branch_name>
```
At this point, GitHub shall have a remote branch that tracks your local branch. You
should be able to see the remote branch by going to the git repo and look for "insights/network"

### 8. Submit a pull request (PR)
This step is done on GitHub site. At this stage, only fill in the title of the PR to
indicate what this branch is for. No content is needed.

From now on, all subsequent commits and pushes to `<branch_name>` will be staged in this
PR and when the PR is merged to the main branch. All commits in this PR will be
squashed. Then, you will write a meaningful title and contents documenting the changes,
use cases and notes of this PR.

## Optional packages
- The [Reference Forward Model](http://eodg.atm.ox.ac.uk/RFM/) (RFM) is provided optionally as
a tool to generate opacity tables. The source code of this package is not publically available.
Please contact the original author or Cheng Li to obtain access. The build process turns on RFM
by default, to turn off building RFM, use
```
cmake -DRFM=OFF ..
```

## Large file storage
This repo uses Git Large File (lfs) Storage to store opacity data.
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

## Quick tips
- undo a "git add"
```
git reset <file>
```
