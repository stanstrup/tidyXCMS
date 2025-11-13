#!/bin/bash
# delete-bad-tags.sh - Delete incorrect version tags from remote
set -e

echo "Deleting incorrect version tags from GitHub..."
echo ""

# Fetch all tags from remote first
echo "Fetching all tags from remote..."
git fetch --tags origin

# List of tags to remove
TAGS_TO_REMOVE="1.0.0 1.0.1 1.0.2 1.0.3 1.0.4 1.0.5 1.0.6 v1.0.0 v2.0.0 v2.0.1 v2.0.2 v2.0.3 2.0.0 2.0.1 2.0.2 2.0.3 semantic-release-1.0.0"

# Find any other semantic-release-* tags
echo "Checking for additional semantic-release-* tags..."
SEMANTIC_TAGS=$(git tag | grep "^semantic-release-" || true)
if [ -n "$SEMANTIC_TAGS" ]; then
  # Add them to the list (avoid duplicates)
  for tag in $SEMANTIC_TAGS; do
    if ! echo "$TAGS_TO_REMOVE" | grep -q "$tag"; then
      TAGS_TO_REMOVE="$TAGS_TO_REMOVE $tag"
    fi
  done
fi

echo "Tags to remove: $TAGS_TO_REMOVE"
echo ""

# Delete tags from remote
echo "Deleting tags from remote repository..."
for tag in $TAGS_TO_REMOVE; do
  if git ls-remote --tags origin | grep -q "refs/tags/$tag$"; then
    echo "Deleting remote tag: $tag"
    git push origin ":refs/tags/$tag" || echo "  Failed to delete $tag (may not exist or no permission)"
  else
    echo "Tag $tag not found on remote"
  fi
done

# Delete local tags
echo ""
echo "Deleting local tags..."
for tag in $TAGS_TO_REMOVE; do
  if git rev-parse "$tag" >/dev/null 2>&1; then
    git tag -d "$tag"
    echo "  ✓ Deleted local tag: $tag"
  fi
done

echo ""
echo "✓ Cleanup complete! Verify at: https://github.com/stanstrup/tidyXCMS/releases"
echo "Run 'git fetch --tags --prune origin' to sync your local tags"
