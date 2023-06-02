# How to use git for this repo 📌

We adopt the idea of **linear history** and a **squash merging** approach in this repository, meaning there is only one permanent branch (`main`), and the only way to push changes to main is by submitting a Pull Request (PR). The main branch is protected to prevent direct pushes. **A linear history ensures that the main branch remains clean and organized**. Squash merging means that **the smallest unit of change is a PR, rather than a commit**.

This workflow differs from some individual workflows where the smallest unit is usually a commit. For collaborative projects, commits can be too fine-grained and don't track issues effectively. Our aim is to ensure that each stage in the history solves a problem that can be traced back, providing context for that problem. In other words, development is **_issue-driven_**. The git workflow recommended for this repository goes as follows:

📍 Step 1. Submit an issue ticket on the GitHub website:

Before starting any work, create an issue ticket on GitHub to describe the problem or task you want to address.

- If you plan to solve the issue yourself, a brief title without extensive details is sufficient.
- If you want someone else to solve it, provide a more detailed explanation.

📍 Step 2. Create a new branch locally

Once you have identified an issue or task, create a new branch locally using the following command:

```bash
git checkout -b <username/issue_description>
```

The branch name should start with your username, followed by a slash `/` and then a brief description of the nature of the issue or task you're working on (for example, `/add_cpptest_case`). This naming convention helps keep the branch names consistent and makes it easier to identify the purpose of each branch, as well as who is working on it.

📍 Step 3. Use the .gitignore file

The `.gitignore` file helps keep your working directory clean. Each folder can have its own `.gitignore` file, which lists files that should not be tracked by the git system. For example, model output files should not be added to git. Ideally, when you run `git status`, there should be no untracked files in your working directory. This state is called a "clean workspace". When you sign off the work for today, please make sure that your workspace is clean. You do not need to worry about staging incomplete work. For your last commit of the day, use a meaningful description message because it is unsure when you will pick it up again.

📍 Step 4. Add changed files to git

If you put all untracked files in `.gitignore`, when you want to pause your work, you can simply use the following command to add your changes to git:

```bash
git add .
```

This command adds all modified files to git (not yet commited), excluding the files listed in the `.gitignore` file.

📍 Step 6. Commit your changes locally

After adding the changed files, use `git status` to review the modifications. If you accidentally added files that you don't want to include, you can undo the add by `git reset <file>`. Then, commit your changes locally with a message:

```bash
git commit -m "<message>"
```

The content of the message is not crucial at this stage since all the commits within a PR will be squashed later into a single message that you'll write later. You can use a meaningful message like "work on XXX" or "working XXX," or a generic message like "wip" (work in progress).

📍 Step 7. Upload your branch to GitHub

The previous command only commits the changes locally. To push your changes to GitHub, use:

```bash
git push origin <username/issue_description>
```

This command pushes your branch to the remote repository on GitHub. You should be able to see the remote branch by visiting the repository's page and looking for "insights/network."

📍 Step 8. Submit a pull request (PR)

This step is performed on the GitHub site. At this stage, only provide a title for the PR to indicate the purpose of the branch. No additional content is required.

All subsequent commits and pushes to `<username/issue_description>` will be included in this PR. When the PR is merged into the `main` branch, all the commits in the PR will be squashed. At that point, you can write a meaningful title and description that document the changes, issues, usage, and any additional notes related to the PR.

📍 Step 9. Update your local branch

After the PR merge, update your local branch by:

```bash
git checkout main
git fetch
git rebase origin/main
git branch -D <username/issue_description>
```

where the `<username/issue_description>` refers to the branch you created in Step 2. This command will delete the local branch. You can then continue to start a new cycle and work on new issues.

---

Following this workflow ensures a systematic approach to contributions, promotes collaboration, and maintains a clean and organized codebase.
