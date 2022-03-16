# Version 0.5.0

This is the first submission of `pmxpartab` to CRAN.

Addressed the following comments from CRAN reviewer:

* From: Julia Haider <julia.haider@wu.ac.at> on 14-Mar-2022:

  - Please add \value to .Rd files regarding exported methods and explain the
    functions results in the documentation. Please write about the structure of
    the output (class) and also what the output means. (If a function does not
    return a value, please document that too, e.g. \value{No return value, called
    for side effects} or similar)
    
    Missing Rd-tags:
    * fpval.Rd: \value
    * pmxpartab.Rd: \value

Also, changed license to MIT.

## Test environments

* Local:
  - Windows 10: R 4.1.3 x86_64-w64-mingw32/x64 (64-bit))
  - Ubuntu Linux 20.04.4 LTS: R 4.1.3 x86_64-pc-linux-gnu (64-bit)
* win-builder:
  - R version 4.1.3 (2022-03-10) x86_64-w64-mingw32 (64-bit)
    - 1 NOTE:
      * New submission
  - R Under development (unstable) (2022-03-12 r81880 ucrt) x86_64-w64-mingw32 (64-bit)
    - 1 NOTE:
      * New submission

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

There are currently no downstream dependencies for this package.

