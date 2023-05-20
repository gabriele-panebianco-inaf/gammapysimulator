#!/usr/bin/env bash
: '
    License     
'

set -e

# Identify directory
SCRIPT_DIRECTORY="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

# Remove __pycache__ directories for a clean install
find $SCRIPT_DIRECTORY -type d -name __pycache__ -exec rm -r {} \+

# Remove package.egg-info directory for a clean install
find $SCRIPT_DIRECTORY -type d -name *.egg-info -exec rm -r {} \+

# Remove .pytest_cache directory for a clean install
find $SCRIPT_DIRECTORY -type d -name .pytest_cache -exec rm -r {} \+

# Install package in development mode
pip install -e .
