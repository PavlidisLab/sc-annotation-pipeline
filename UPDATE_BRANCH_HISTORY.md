git checkout development
git reset --hard origin/development

# Updating Your Local Repository After Branch History Reset

The `main` and `development` branches have been force-updated to a new history. To continue working with the latest code, all users must update their local repositories as follows:

## 1. Fetch and Prune Remotes
Fetch all updates and remove any stale remote-tracking branches:

```
git fetch origin --prune --force
```

## 2. Verify Remote Branch History
Check that `origin/main` and `origin/development` show the new commit history:

```
git log origin/main --oneline -n 5
git log origin/development --oneline -n 5
```

If these do not show the new history, try removing and re-adding the remote:

```
git remote remove origin
git remote add origin <REPO_URL>
git fetch origin --prune --force
```

## 3. Reset Local Branches to Match Remote
**Warning:** This will overwrite your local `main` and `development` branches. Backup any local changes first!

```
git checkout main
git reset --hard origin/main

git checkout development
git reset --hard origin/development
```

## 5. Notes
- If you have local changes you want to keep, commit or stash them before running `reset --hard`.
- If you have forks or feature branches based on the old history, you may need to rebase or recreate them.
- If you encounter errors, contact the repository maintainer for help.

---

**Why is this necessary?**

The commit history of `main` and `development` was force-updated to resolve divergence and make `my-salvaged-changes` the new base. All users must sync to this new history to avoid conflicts.
