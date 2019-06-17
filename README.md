# <img src="https://twiki.cern.ch/twiki/pub/BEABP/Logos/OMC_logo.png" height="28"> Beta-Beat Source

This is the python-tool package of the optics measurements and corrections group (OMC).

If you are not part of that group, you will most likely have no use for the codes provided here, 
unless you have a 9km wide accelerator at home.
Feel free to use them anyway, if you wish!

## Documentation

- Autgenerated docs via ``sphinx`` can be found on <https://pylhc.github.io/Beta-Beat.src/>.
- General documentation of the OMC-Teams software on <https://twiki.cern.ch/twiki/bin/view/BEABP/OMC>

## Getting Started

### Prerequisites

The codes use a multitude of packages as can be found in the [requirements.txt](requirements.txt).

Important ones are: ``numpy``, ``pandas`` and ``scipy``.

### Installing

This package is not deployed, hence you need to use the standard git-commands to get a local copy.

## Description

This is the old repository ([new one](https://github.com/pylhc/omc3)) of the codes,
written for python 2.7.  


## Quality checks

### Tests

The following tests are run automatically after each commit via 
[Travis-CI](https://travis-ci.org/pylhc/Beta-Beat.src):

- Pytest unit tests
- Accuracy tests
- Regression tests

### Maintainability

- Additional checks for code-complexity, design-rules, test-coverage, duplication on 
[CodeClimate](https://codeclimate.com/github/pylhc/Beta-Beat.src)


## Authors

* **pyLHC/OMC-Team** - *Working Group* - [pyLHC](https://github.com/orgs/pylhc/teams/omc-team)

<!--
## License
This project is licensed under the  License - see the [LICENSE.md](LICENSE.md) file for details
-->