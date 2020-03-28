# Change Log
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/)
and this project (attempts to) adhere to [Semantic Versioning](http://semver.org/).

## [0.5.2] - 2020-03-28
- Remove redundencies in b6 when loading (keep the one with the lowest e-value)

## [0.5.1] - 2020-03-25
- Add "SkipGenes"

## [0.5.0] - 2020-03-25
- Add an "--update" option to tax_collector
- Update the README a tiny bit
- Change default to group rather than species

## [0.4.3] - 2020-03-09
- so dumb

## [0.4.2] - 2020-03-09
- change how functional_annotation.py works

## [0.4.1] - 2020-03-09
- fix dumb bug with wrong name of functional_annotation.py

## [0.4.0] - 2020-03-09
- tax_collector now correctly produces "tsv" files
- added "functional_tax.py"; we'll see how it does

## [0.3.0] - 2020-03-05
- Allow tax_collector to accept an .stb file

## [0.2.0] - 2020-03-04
- Changes to the README to make it actually reflect what it can do
- Updates to allow it to process files created using the diamond pipeline

## [0.1.0] - 2017-06-07
### Fixed
- program "make_Tdb.py" works again (bug fix with regards to pandas merge update)
- TaxIDs are now internally stored as floats in some places (when np.nan values can be there)
- very minor documentation added

## [0.0.3] - 2017-11-29
### Added
- program "make_Tdb.py" works again

## [0.0.2] - 2017-11-27
### Added
- program "tax_collector.py" added to bin

## [0.0.1] - 2017-11-07
### Added
- the beast is born
