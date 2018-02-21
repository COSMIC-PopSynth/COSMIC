#!/bin/bash

# enable strict test flags
if [ "$STRICT" = true ]; then
    _strict="-x --strict"
else
    _strict=""
fi

coverage run --append `which runFixedPop` --help
