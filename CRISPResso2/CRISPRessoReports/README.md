# CRISPRessoReports
Git Subtree Standard Operating Procedure (SOP)

Created: February 3, 2023

Introduction
In a few of our repositories we use git subtrees to share common files shared across multiple repositories. The goal of using subtrees is to reduce duplicated code.

Process
If you have a repository where you would like to add a subtree, run the following command from the parent repo (wherein you would like the external child repo to be stored).

git subtree add --prefix {path to where subtree will be stored} --squash {url of child repo} {branch of child repo}

Some notes about the above command:
- You can’t add a subtree to an empty repo.
- If you don’t specify the --squash parameter then all of the commits from the child repo will be brought into the parent repo history.
- When you do specify the --squash parameter then there will be two commits in the parent repo history.
- If you have the url of the child repo added as a remote, you can use that instead.
- You should now see the files in the child repo at the path you specified in the parent repo.

If you edit the files in the child repo (from within the parent repo), you can push the changes you have made up to the child repo with this command.

git subtree push --prefix {path to where subtree is stored} --rejoin {url of child repo} {branch of child repo}

Some notes about the above command:
- You should separate your commits such that there are not files staged in both child and parent.
- If there are commits in the child repo that you need in the parent repo, you can use this command.

git subtree pull --prefix {path to where subtree is stored} --squash {url of child repo} {branch of child repo}

Some notes about the above command:
- A merge conflict can happen here, so be warned!

If you are doing this regularly, the commands above can be cumbersome. You can create an alias in git to make these commands easier! Edit the .git/config file in your repo and add the following lines to it.

[alias]

The acronym stands for "CrispressoReports Add"

cra = "!f() { git subtree add --prefix CRISPRessoWEB/CRISPRessoReports https://github.com/edilytics/CRISPRessoReports.git $1 --squash; }; f"

The acronym stands for "CrispressoReports Update"

cru = "!f() { git subtree pull --prefix CRISPRessoWEB/CRISPRessoReports https://github.com/edilytics/CRISPRessoReports.git $1 --squash; }; f"

The acronym stands for "CrispressoReports Push"

crp = "!f() { git subtree push --prefix CRISPRessoWEB/CRISPRessoReports https://github.com/edilytics/CRISPRessoReports.git $1; }; f"


