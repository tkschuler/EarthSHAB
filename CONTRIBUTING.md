# Contributing to EarthSHAB

Thanks for contributing! This guide covers how to set up the project, make a
change, and open a pull request.

## Workflow

1. **For non-trivial changes, open (or comment on) an issue first** describing
   the approach. A quick design comment before coding avoids wasted work and
   gives a thread to discuss in. Reference it from your PR with `Closes #<n>`.
2. **Fork the repo** and clone your fork.
3. **Branch off `main`**, leading with the issue number you're addressing:
   `issue-<#>-<short-description>` (e.g. `issue-5-spatio-temporal-interp`); the
   shorter `<#>-<short-description>` is fine too, and is what GitHub's
   *Create a branch* button generates from an issue. For a trivial change with no 
   issue, a short descriptive name (e.g. `fix-readme-typo`) is fine.
4. **Make your change** with tests and docs (see below).
5. **Open a PR**, with `main` as the base branch. Push more commits
   to the same branch as you iterate — the PR updates automatically.
6. Mark it **Ready for review** once CI is green.

Branch always off the latest `main`:

```bash
git checkout main && git pull
git checkout -b issue-5-my-change         # issue-<#>-<short-description>
# ...edit, commit...
git push -u origin issue-5-my-change
```

Open the PR (base `main`):

```bash
gh pr create --base main --draft \
  --title "Short description" --body "Closes #<n>"
```

When opening your PR, tick **"Allow edits by maintainers"** so small
fixes can be pushed to your branch during review.

Keep your branch current if `main` moves while your PR is open:

```bash
git fetch origin && git rebase origin/main      # or: git merge origin/main
```


## Tests

```bash
pytest tests/
```

Add tests for new behavior. Pure, deterministic logic should have unit tests;
mirror the existing patterns (e.g. `tests/test_wind_interpolation.py`, and the
synthetic forecast fixtures in `tests/tools/make_dummy_forecasts.py`).

## Documentation & changelog

- Update the relevant docs under `docs/source/` when you change behavior or
  config. You can build them locally:
  `sphinx-build -W -b html docs/source docs/html`.
- Add a short bullet to the `Changelog`.
