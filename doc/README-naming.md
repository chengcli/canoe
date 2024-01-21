
## Naming Conventions

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

### How to name folders and files
- use a singular simple noun for folder names
- avoid compound nouns for folder names
- you can use compound nouns or phrases for files names
- individual words in a compound nouns for phrase should be concatenated by underscore

### How to name variables and classes
- use low case letters for variables
- you can use compound nouns or phrases
- compound nouns for phrase should be concatenated by underscore. This is called the *snake case*
- if a compound noun is used for class name, capitalize each word in the compound noun.
  This is called the *upper camel case*
- if a variable is a private member of a class, append the variable name with an
  underscore

### How to name funcions
- functions names are usualy verbal phrases
- standalone functions should use *snake case*
- public class member functions should use *upper camel case*
- private class member functions should use *lower camel case* in which the first word
  in a phase is not capitalized and the rest words are capitalized
