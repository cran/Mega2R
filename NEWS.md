# Version 1.0.5 (2019-02-27)

-- Bug fix: The 'mkfam' function now correctly extracts the case/control trait when the database contains more than one phenotype.  It worked correctly previously in our simple example data set where there was only one phenotype in the Mega2 database.

-- Mega2pedgene adjusted to allow specification of the trait name.

# Version 1.0.4 (2018-06-18)

-- Improvements to use compressed data created by Mega2 version 5.0.0 or higher.

# Version 1.0.3 (2018-05-22)

-- Removed strict dependency on GenABEL because GenABEL has been archived.

# Version 1.0.2 (2018-04-03)

- Bug fix: The init_pedgene function now sets up the trait and pedigree structure correctly.

# Version 1.0.0 (2017-08-22)

- Initial CRAN release
