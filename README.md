# Microwaves

This package is ultimately intended to provide a library of functions useful
for RF and microwave engineering. Currently, the main feature being
implemented is a touchstone file reader and plotting recipes for the read
data.

There is currently a semi-functional reader named `read_touchstone()`,
however this does not work for all touchstone files, and it will be
superseded soon by interfacing with FileIO so users can simply call `load()`.

Additionally, functionality for network manipulation and basic network
synthesis is planned. Any other useful ideas/features people can come up
with are more than welcome as well.

If it is not apparent, this package is in the early stages of its development
and its API should not be considered stable.