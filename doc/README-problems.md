
# Common problems

## Cannot fetch RFM
- Ask chengcli@umich.edu for access to RFM (privately hosted)
- Give **full access** to your Github access token (classic)
- Set environment variable `GH_ACCOUNT` to your Github username
- Set environment variable `GH_TOKEN` to your Github access token (classic)
- Source your `.bashrc` or `.zshrc` file

## Cannot fetch RFM in Github actions
- Set Github action secrets
    1. Go your repo -> `settings` -> `Secrets and variables` -> `Actions`.
    1. Set a `New repository secret` with name `ACCESS_TOKEN`.
    1. The value of the secret is your github access toekn.
