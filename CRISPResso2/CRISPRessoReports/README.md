# CRISPRessoReports

This repo holds the shared reports code and HTML templates that are used by the CLI and web projects.

Take care when committing into these files as not to mix unrelated git histories.

## How do I work with this repo?

Step 1 only needs to be done once per cloned repo, the other steps will need to be dne more frequently.

1. Add the remote to this repo to the "parent" repo (i.e. `CRISPResso2` or `C2Web`). **You should only have to do this once.**

``` shell
git remote add reports https://github.com/edilytics/CRISPRessoReports.git
```

*Note:* you can change the `reports` above to whatever you would like, just don't forget to replace `reports` with it in the other steps.

2. In the parent repo, fetch the latest changes from `CRISPRessoReports`.

``` shell
git fetch reports
```

3. Checkout a `CRISPRessoReports` branch.

``` shell
git checkout -b master-reports reports/master
```

*Note:* you can obviously change the names of these branches (or select a branch other than `master`). I would recommend naming all branches from `CRISPRessoReports` with the same prefix or suffix (i.e. `*-reports`) because it can get confusing to keep all of the branches straight...

4. Read `CRISPRessoReports` into the parent repo as a subtree. This is when the code in `CRISPRessoReports` will finally be added to the parent repo.

**Very important**, switch to the branch in the parent repo where you want the code to be added! For example:

``` shell
git checkout <feature-branch>
```

Then, you will read the commits into wherever `CRISPRessoReports` are stored for that repo. **Note:** you should only have to do this if `CRISPRessoReports` has not been added to the parent repo, if it is already there, do not repeat this step.

``` shell
git read-tree --prefix=CRISPResso2/CRISPRessoReports -u master-reports
```

Run `git status` and you should see the files added!

5. Stage and commit your files as you normally would.

### How do I pull commits that are in `CRISPRessoReports` into the parent repo?

1. In the parent repo, switch to the `<feature-branch>-reports` branch.

``` shell
git checkout <feature-branch>-reports
```

2. Pull the changes from `CRISPRessoReports`.

``` shell
git pull
```

You should see the updates that you are looking for.

3. **Very important**, switch back to whichever branch you are working on in the parent repo.

``` shell
git checkout <feature-branch>
```

4. Merge the changes in and resolve any merge conflicts.

``` shell
git merge --squash -Xsubtree="CRISPResso2/CRISPRessoReports" --no-commit --allow-unrelated-histories <feature-branch>-reports
```

*Note:* You may need to change the value of the `-Xsubtree` parameter to match where `CRISPRessoReports` is located in the parent repo.

5. Commit your changes and resolve merge conflicts.

Also, note that the default commit message may have a summary of all commits, please shorten it to be descriptive of the most recent changes.

### How do I push commits that are in my parent repo's `CRISPRessoReports` into the shared `CRISPRessoReports` repo?

1. In the parent repo, switch to (or create) the branch on `CRISPRessoReports` that will have the changes you push.

If you are creating a new branch based off of `CRISPRessoReports` master, run this to switch to the reports master branch:

``` shell
git checkout reports/master
```

Then, run to actually create (and switch to) the branch that you will be working with:

``` shell
git checkout -b <feature-branch>-reports
```

Or if you would like to push to an existing branch on `CRISPRessoReports`, run this:

``` shell
git checkout <feature-branch>-reports
```

2. Merge the changes in and resove any merge conflicts.

``` shell
git merge --squash -Xsubtree="CRISPResso2/CRISPRessoReports" --no-commit --allow-unrelated-histories <feature-branch>
```

*Note:* `<feature-branch>` is the branch of the parent repo that contains the changes inside the `CRISPRessoReports` sub-directory.

3. Push to `CRISPRessoReports`.

``` shell
git push
```

4. Switch back to your branch on `CRISPResso` or `C2Web`.

``` shell
git checkout <feature-branch>
```

### I am working on a feature that requires changing `CRISPRessoReports`, what do I do?

If a feature that you are working on requires changes to CRISPRessoReports, you will need to perform a few steps to get setup.

1. Create a feature branch in the parent repo, based on the parent repo master.

``` shell
git checkout -b <feature-branch>
```

2. Create a feature branch on `CRISPRessoReports`.

Checkout your local `CRISPRessoReports` master branch.

``` shell
git checkout master-reports
```

Pull in the latest changes.

``` shell
git pull
```

Create the `CRISPRessoReports` feature branch based on `reports/master`.

``` shell
git checkout -b <feature-branch>-reports reports/master
```

*Note:* if your branch is named `cool-feature` in the parent repo, then follow the convention of naming the corresponding `CRISPRessoReports` branch `cool-feature-reports`.

If you run `git status` at this point you should see any directories in the parent repo as untracked files, this is normal and expected.

3. Switch back to the feature-branch in the parent repo, and develop your changes.

``` shell
git checkout <feature-branch>
```

*Note:* you can mingle your changes in `CRISPRessoReports` and the parent repo in the same commits.

4. Merge and push your changes up to `CRISPRessoReports`.

Switch to the `<feature-branch>-reports` branch.

``` shell
git checkout <feature-branch>-reports
```

Merge the changes from the parent repo into the `<feature-branch>-reports` branch.

``` shell
git merge --squash -Xsubtree="CRISPResso2/CRISPRessoReports" --no-commit --allow-unrelated-histories <feature-branch>
```

# FAQ

## There are lots of merge conflicts, how do I just accept all of the incoming changes?

If you want to blindly accept all of the incoming changes, you can add the parameter `-Xtheirs` to the `git merge...` command and anything that was a merge conflict before, should now be overwritten by the incoming change.

## I tired of typing `git merge --squash ...`, what can I do?!

Typing out the `git merge...` command is long and is a big pain. Here are some shortcuts to add to your `.git/config` file in order to make this easier to type.

``` git-config
[alias]
    # merge in branch and resolve merge conflicts
    m = "!f() { git merge --squash -Xsubtree='CRISPResso2/CRISPRessoReports' --no-commit --allow-unrelated-histories $1; }; f"
    # merge in branch and accept all of the incoming changes
    mt = "!f() { git merge --squash -Xtheirs -Xsubtree='CRISPResso2/CRISPRessoReports' --no-commit --allow-unrelated-histories $1; }; f"
```

Now you can just run `git m <feature-branch>` to merge `<feature-branch>` into your current branch. Or run `git mt <feature-branch>` to accept all of the incoming changes.

# Sources and helpful links

- This method was heavily based off of what was described in [this blog post](http://johnatten.com/2013/03/16/git-subtree-merge-the-quick-version/)
- If you want to know more about git merging, the manual is [very helpful](https://git-scm.com/book/en/v2/Git-Tools-Advanced-Merging)
- If you messed up (merging the wrong branch), you can undo it using `git reset --hard <branch>`. **Beware:** this can cause you to lose work, so use with care. [Learn more here](https://stackoverflow.com/a/8888015/1947159).
- If you need to rewrite git history, try using [git-filter-repo](https://github.com/newren/git-filter-repo)
- After rewriting git history (from a mirror repo), if you can't push to GitHub, [try this](https://stackoverflow.com/a/34266401/1947159)
