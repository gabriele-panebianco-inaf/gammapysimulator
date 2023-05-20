#######################################################
#
# Authors:
#
#   Gabriele Panebianco <gabriele.panebianco@inaf.it>
#
#######################################################

#!/bin/bash

set -e 

SCRIPT_DIRECTORY="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
export GAMMAPYSIMULATOR=$SCRIPT_DIRECTORY

# Make artifacts directories for test reports
DIRECTORY_JUNIT="$GAMMAPYSIMULATOR/.artifacts/junit-reports"
DIRECTORY_COVERAGE="$GAMMAPYSIMULATOR/.artifacts/coverage-reports"

mkdir -p "$DIRECTORY_JUNIT" "$DIRECTORY_COVERAGE"
rm -f "$DIRECTORY_JUNIT/*" "$DIRECTORY_COVERAGE/*"

cd $SCRIPT_DIRECTORY/test

pytest \
 --junitxml=$DIRECTORY_JUNIT/tests_report.xml \
 --cov=$GAMMAPYSIMULATOR/gammapysimulator \
 --cov-config=.coveragerc \
 -vv \
 --exitfirst \
 -s

coverage xml -o "$DIRECTORY_COVERAGE/coverage_report_xml/coverage_report.xml"
coverage html -d "$DIRECTORY_COVERAGE/coverage_report_html"

echo "Coverage report: $DIRECTORY_COVERAGE"
echo "Test results: $DIRECTORY_JUNIT"

