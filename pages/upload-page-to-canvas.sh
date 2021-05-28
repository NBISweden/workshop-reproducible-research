#!/bin/bash

# Script that takes a markdown file as input and creates or updates a page for
# the course based on its filename, as applicable.

# The markdown document to be uploaded should be supplied as the first argument
# to the script. The title of the upload can be changed in Canvas after the
# first upload to something more suitable, which won't affect subsequent
# updates.

# Images are linked to their specified location on GitHub (the `GITHUB`
# variable) and should be referenced to as e.g. `![](images/<some-image.png>)`
# in the markdown.
#
# Links to pages on Canvas can be given using the URL to the page, but replacing
# the course ID with the string 'COURSE_ID', which will then automatically build
# the correct link by using the $COURSE_ID variable defined below.

# Input parameters
MARKDOWN=$1
COURSE_ID=51980

# General parameters
TOKEN="$HOME/.canvas-api-token"
API="https://uppsala.instructure.com/api/v1/courses"
PAGE=$(basename $MARKDOWN | sed 's/.md//g')
HTML=$(basename $MARKDOWN | sed 's/.md/.html/g')

# Current Git branch and course image path
BRANCH=$(git branch | sed -n -e 's/^\* \(.*\)/\1/p')
GITHUB="https:\/\/raw\.githubusercontent\.com\/NBISweden\/workshop-reproducible-research\/$BRANCH\/pages\/"

# Convert using Pandoc
docker run --rm \
    --volume "`pwd`:/data" \
    --user `id -u`:`id -g` \
    pandoc/latex $MARKDOWN --output="$HTML"

# Add images from GitHub and course ID for links
cat "$HTML" \
    | sed "s/\(src=\"\)\(images\/\)/\1$GITHUB\2/g" \
    | sed "s/COURSE_ID/$COURSE_ID/g" \
    > tmp.html
mv tmp.html "$HTML"

# Create or update page curl PUT
echo "Uploading \`$HTML\` ..."
curl -X PUT \
    "$API/$COURSE_ID/pages/$PAGE" \
    --header "Authorization: Bearer $(cat ~/.canvas-api-token)" \
    --data-urlencode wiki_page[body]="$(cat $HTML)" \
    --silent --show-error \
    > /dev/null

# Delete rendered HTML
rm $HTML
