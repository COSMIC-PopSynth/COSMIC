#!/bin/bash

# enable strict test flags
if [ "$STRICT" = true ]; then
    _strict="-x --strict"
else
    _strict=""
fi

coverage run --append `which cosmic-pop` --help
coverage run -m py.test -v -r s ${_strict} cosmic/
