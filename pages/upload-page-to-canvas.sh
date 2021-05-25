#!/bin/bash

# Script that takes a markdown file as input and creates or updates a page for
# the course based on its filename, as applicable. The markdown document to be
# uploaded should be supplied as the first argument to the script. The title of
# the upload can be changed in Canvas after the first upload to something more
# suitable, which won't affect subsequent updates.

# Parameters
MARKDOWN=$1
COURSE_ID=51980
TOKEN="$HOME/.canvas-api-token"
API="https://uppsala.instructure.com/api/v1/courses"
PAGE=$(basename $MARKDOWN | sed 's/.md//g')
HTML=$(basename $MARKDOWN | sed 's/.md/.html/g')
PAGE_WIDTH=4

# Convert using Pandoc
docker run --rm \
    --volume "`pwd`:/data" \
    --user `id -u`:`id -g` \
    pandoc/latex $MARKDOWN -o $HTML

# Add maximum page width
CONTENT="$(echo \<div class="col-lg-$PAGE_WIDTH"\>; cat $HTML; echo \</div\>)" \
    > /dev/null

# Check if current page already exists
curl -X GET \
    "$API/$COURSE_ID/pages/$PAGE" \
    --header "Authorization: Bearer $(cat $TOKEN)" \
    --silent --show-error \
    > /dev/null

# Create or update page, as applicable
if [ $? -eq 0 ]; then
    echo "Page \`$PAGE\` already exists: updating it"
    curl -X PUT \
        "$API/$COURSE_ID/pages/$PAGE" \
        --header "Authorization: Bearer $(cat ~/.canvas-api-token)" \
        --data wiki_page[body]="$(echo $CONTENT)" \
        --silent --show-error \
        > /dev/null
else
    echo "Page \`$PAGE\` does not exist: creating it"
    curl -X POST \
        "$API/$COURSE_ID/pages/$PAGE" \
        --header "Authorization: Bearer $(cat ~/.canvas-api-token)" \
        --data wiki_page[title]="$PAGE" \
        --data wiki_page[body]="$(echo $CONTENT)" \
        --silent --show-error \
        > /dev/null
fi

# Delete rendered HTML
rm $HTML
