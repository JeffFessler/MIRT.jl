To make edits or additions to `MIRT.jl`
first clone `MIRT.jl`
or use `]dev MIRT` to make a local version.

Then it probably is best to make a new branch (called `try`) here
for trying out changes.

* `git checkout -b try`
* `git push origin try`
* `git branch -a`

Add files as needed.
Be sure to comment out "using MIRT"
and update `z-test.jl` and `z-list.jl`

See https://github.com/JeffFessler/MIRT.jl/pull/15/files
for an example of an added file and the corresponding changes.

* `git add .`
* `git commit`
* `git push --set-upstream origin try`

Wait for build to pass,
then check that code coverage is 100%.
Edit as needed.
When ready, submit pull request on github.com.

After pull request is approved:
* `git pull`
* `git checkout master
* `git pull
