# Writing with Rmarkdown

We can create versioned (ie in Github), live (ie output directly from data with R), documents in various formats (ie pdf, html, docx) with formatting, figures, tables, equations and references. All using [Rmarkdown](http://rmarkdown.rstudio.com). This is in the vein of generating truly reproducible research.

## Dependencies using renv
[renv](https://github.com/rstudio/renv) is a package manager maintained by the RStudio folks.

- To depend on a new package, just add `require(whatever)` statements to your R code, and then run `renv::snapshot()`.
  This ensures the `renv.lock` file contains the exact dependency version you downloaded version.
- Things should Just Work if you open the project in RStudio, but if not you can always run `renv::restore()`. If `renv`
  isn't available, just do the usual `install.packages('renv')`

## Prove your script works using CircleCI

Every time you push a git commit, CircleCI will notice and kick off a 'build' job.  It will clone this repository, as if
it was a competely new team-mate starting from scratch, install all the dependencies and try to run the script (as
defined in [`.circleci/config.yml`](.circleci/config.yml)). Assuming your code is error-free it'll post back a nice
green tick, otherwise you get a red X.

Check out the [CircleCI 'Artifacts' tab](https://app.circleci.com/pipelines/github/iamdanfox/rmarkdown-example/21/workflows/fb58c7a9-0cc3-4c57-8c5d-494c3b45ea92/jobs/22/artifacts) to see Rmarkdown rendered to HTML.

