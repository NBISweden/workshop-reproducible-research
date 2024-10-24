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

# Links to pages on Canvas can be given using the URL to the page, but replacing
# the course ID with the string 'COURSE_ID', which will then automatically build
# the correct link by using the $COURSE_ID variable defined below. The COURSE_ID
# needs to be updated between course rounds, which is done by editing the
# `pages/.course_id` file. You can also just provide the link using the markdown
# format: `[link text](<page-name>)`.

# This script is run automatically in a GitHub Actions workflow which you can
# find at `.github/workflows/canvas-upload.yml`. There is only one page that it
# will run on, and that is the main landing page used as an example above, while
# the rest of the course materials lives on GitHub pages.

# Input parameters
MARKDOWN=$1

if [ "$2" == "" ]; then
  TOKEN=$(cat "$HOME/.canvas-api-token")
else
  TOKEN="$2"
fi

# General parameters
API="https://uppsala.instructure.com/api/v1/courses"
PAGE=$(basename $MARKDOWN | sed 's|.md||g')
HTML=$(basename $MARKDOWN | sed 's|.md|.html|g')

# Get current branch and build GitHub address
BRANCH=$(git branch | sed -n -e 's|^\* \(.*\)|\1|p')
GITHUB="https://raw.githubusercontent.com/NBISweden/workshop-reproducible-research/$BRANCH/pages/"

# Get the appropriate course ID from the current branch
# (set to development ID on feature branches for testing)
COURSE_ID=$(grep $BRANCH pages/.course_id | cut -f2 -d ':')
if [ "$COURSE_ID" == "" ]; then
    COURSE_ID=54324
fi

# Convert using Pandoc
echo "Rendering \`$MARKDOWN\` ..."
docker run --rm \
    --volume "`pwd`:/data" \
    --user `id -u`:`id -g` \
    pandoc/latex $MARKDOWN --mathjax --output="$HTML"

# Add images from GitHub and course ID for links
cat "$HTML" \
    | sed "s|\(src=\"\)\(images/\)|\1$GITHUB\2|g" \
    | sed "s|COURSE_ID|$COURSE_ID|g" \
    | sed "s|GITHUB_BRANCH|$BRANCH|g" \
    > tmp.html
echo '<div class="content" style="max-width:800px">' \
    | cat - tmp.html \
    > tmp2.html
echo "</div>" >> tmp2.html
mv tmp2.html "$HTML"
rm tmp.html

# Create or update page curl PUT
echo "Uploading \`$HTML\` to course ID $COURSE_ID ..."
curl -X PUT \
    "$API/$COURSE_ID/pages/$PAGE" \
    --header "Authorization: Bearer $TOKEN" \
    --data-urlencode wiki_page[body]="$(cat $HTML)" \
    --silent \
    --show-error \
    --fail \
    --output /dev/null

# Store curl exit code
CURL_EXIT_CODE=$?

# Remove rendered HTML page
rm $HTML

# Exit with curl exit code
if [ $CURL_EXIT_CODE != 0 ]; then
    exit $CURL_EXIT_CODE
fi
