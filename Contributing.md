# Contributing guide

## Developer conventions

Developers are expected to follow the conventions described here.

### Issues
Should follow the bug/refactor/feature template.

### Branches
Development should happen primarily on a separate branch and merged into `ARCCollaboration` (the main development branch) via a pull request (PR) following the PR template. Changes should not be pushed directly to the main development branch (this cannot be enforced in private repositories).

Branches should be named as `developername-issue_number_(s)-some_description`. E.g. `wgraham-19_20-improve_build_system` would be a good name for a branch addressing issues 19 and 20 which entail improvements to the build system.

### Pull requests
At least the approval of one other developer is required before merging a PR. Should the reviewer have no comments or changes they'd like to request, they can merge immediately after approval. Alternatively, they can approve and make some (non-critical) comments for the developer, who can merge after considering these. PRs should not be merged if the reviewer commented or requested changes only.

Development branches should be deleted after merging. Avoid squashing the commits on merging, create a merge commit instead.

* **C++ style guide**  
We loosely aim to move the code base towards adhering to [Google's C++ style guide](https://google.github.io/styleguide/cppguide.html) in order to increase legibility (at least newly introduced functionality should follow this.)
