#!/bin/sh
# here added ross march 20
HERE=`dirname $0`
echo "We are $HERE"
cd $HERE
export PYTHON_EGG_CACHE=$HERE
python ./scripts/check_python.py
[ $? -ne 0 ] && exit 1

SAMPLES="
    external_service_types_conf.xml.sample
    datatypes_conf.xml.sample
    reports_wsgi.ini.sample
    tool_conf.xml.sample
    tool_data_table_conf.xml.sample
    universe_wsgi.ini.sample
    tool-data/shared/ucsc/builds.txt.sample
    tool-data/*.sample
    static/welcome.html.sample
"

# Create any missing config/location files
for sample in $SAMPLES; do
    file=`echo $sample | sed -e 's/\.sample$//'`
    if [ ! -f "$file" -a -f "$sample" ]; then
        echo "Initializing $file from `basename $sample`"
        cp $sample $file
    fi
done

# explicitly attempt to fetch eggs before running
FETCH_EGGS=1
for arg in "$@"; do
    [ "$arg" = "--stop-daemon" ] && FETCH_EGGS=0; break
done
if [ $FETCH_EGGS -eq 1 ]; then
    python ./scripts/check_eggs.py quiet
    if [ $? -ne 0 ]; then
        echo "Some eggs are out of date, attempting to fetch..."
        python ./scripts/fetch_eggs.py
        if [ $? -eq 0 ]; then
            echo "Fetch successful."
        else
            echo "Fetch failed."
            exit 1
        fi
    fi
fi

python ./scripts/paster.py serve universe_wsgi.ini $@
