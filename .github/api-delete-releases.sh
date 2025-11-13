#!/bin/bash
# api-delete-releases.sh - Delete releases and tags using GitHub API
set -e

# Check for GITHUB_TOKEN
if [ -z "$GITHUB_TOKEN" ]; then
  echo "ERROR: GITHUB_TOKEN environment variable is not set"
  echo "Please set it with: export GITHUB_TOKEN=your_token_here"
  exit 1
fi

REPO="stanstrup/tidyXCMS"
API_URL="https://api.github.com"

echo "Deleting incorrect releases and tags from $REPO..."
echo ""

# List of tags to remove
TAGS_TO_REMOVE="1.0.0 1.0.1 1.0.2 1.0.3 1.0.4 1.0.5 1.0.6 v1.0.0 v2.0.0 v2.0.1 v2.0.2 v2.0.3 2.0.0 2.0.1 2.0.2 2.0.3 semantic-release-1.0.0"

for tag in $TAGS_TO_REMOVE; do
  echo "Processing: $tag"

  # Get release ID for this tag
  RELEASE_ID=$(curl -s -H "Authorization: token $GITHUB_TOKEN" \
    "$API_URL/repos/$REPO/releases/tags/$tag" | \
    grep '"id":' | head -1 | sed 's/.*"id": \([0-9]*\).*/\1/')

  if [ -n "$RELEASE_ID" ] && [ "$RELEASE_ID" != "null" ]; then
    echo "  Found release ID: $RELEASE_ID"

    # Delete the release
    HTTP_CODE=$(curl -s -o /dev/null -w "%{http_code}" \
      -X DELETE \
      -H "Authorization: token $GITHUB_TOKEN" \
      "$API_URL/repos/$REPO/releases/$RELEASE_ID")

    if [ "$HTTP_CODE" = "204" ]; then
      echo "  ✓ Deleted release"
    else
      echo "  ✗ Failed to delete release (HTTP $HTTP_CODE)"
    fi
  else
    echo "  - No release found"
  fi

  # Delete the tag (whether or not there was a release)
  HTTP_CODE=$(curl -s -o /dev/null -w "%{http_code}" \
    -X DELETE \
    -H "Authorization: token $GITHUB_TOKEN" \
    "$API_URL/repos/$REPO/git/refs/tags/$tag")

  if [ "$HTTP_CODE" = "204" ]; then
    echo "  ✓ Deleted tag"
  elif [ "$HTTP_CODE" = "422" ] || [ "$HTTP_CODE" = "404" ]; then
    echo "  - Tag not found"
  else
    echo "  ✗ Failed to delete tag (HTTP $HTTP_CODE)"
  fi

  echo ""
done

echo "✓ Cleanup complete!"
echo "Run 'git fetch --tags --prune origin' to sync your local repository"
