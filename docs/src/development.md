# Development

This page contains some helpful information for the development of TrixiBottomTopography.jl.
Further information about useful tools for package development in Julia can be found on the
[development page](https://trixi-framework.github.io/Trixi.jl/stable/development/) of the Trixi.jl docs.

## Releasing a new version of TrixiBottomTopography

- Check whether everything is okay, tests pass etc.
- Set the new version number in `Project.toml` according to the Julian version of semver.
  Commit and push.
- Comment `@JuliaRegistrator register` on the commit setting the version number.
- `JuliaRegistrator` will create a PR with the new version in the General registry.
  Wait for it to be merged.
- Increment the version number in `Project.toml` again with suffix `-DEV`. For example,
  if you have released version `v0.2.0`, use `v0.2.1-DEV` as new version number.



## Preview the documentation

You can build the documentation of TrixiBottomTopography.jl locally by running
```bash
julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
julia --project=docs --color=yes docs/make.jl
```
from the TrixiBottomTopography.jl main directory. Then, you can look at the html files generated in
`docs/build`.
For PRs triggered from branches inside the TrixiBottomTopography.jl main repository previews of
the new documentation are generated at
`https://trixi-framework.github.io/TrixiBottomTopography.jl/previews/PRXXX`,
where `XXX` is the number of the PR.
Note, this does not work for PRs from forks for security reasons (since anyone could otherwise push
arbitrary stuff, including malicious code).
