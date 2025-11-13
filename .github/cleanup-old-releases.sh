#!/bin/bash
# cleanup-old-releases.sh - Remove incorrect 1.0.x releases and tags
# This script should be run once to clean up the repository
set -e

echo "Cleaning up incorrect 1.0.x releases and tags..."

# List of tags to remove
TAGS_TO_REMOVE="1.0.0 1.0.1 1.0.2 1.0.3 1.0.4 1.0.5 1.0.6 v1.0.0"

# Delete GitHub releases (this also removes the tags)
for tag in $TAGS_TO_REMOVE; do
  echo "Deleting GitHub release and tag: $tag"
  if gh release view "$tag" &>/dev/null; then
    gh release delete "$tag" --yes
    echo "  ✓ Deleted release $tag"
  else
    echo "  - Release $tag not found (already deleted)"
  fi

  # Also delete the tag if it still exists on remote
  if git ls-remote --tags origin | grep -q "refs/tags/$tag$"; then
    git push origin --delete "$tag" || echo "  - Could not delete remote tag $tag"
  fi
done

# Delete local tags if they still exist
for tag in $TAGS_TO_REMOVE; do
  git tag -d "$tag" 2>/dev/null || true
done

echo ""
echo "Cleanup complete! Only 0.99.x releases should remain."
echo "Verify at: https://github.com/stanstrup/tidyXCMS/releases"
