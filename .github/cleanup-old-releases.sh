#!/bin/bash
# cleanup-old-releases.sh - Remove incorrect 1.0.x releases and tags
# This script should be run once to clean up the repository
# Requires: gh CLI (GitHub CLI) to be installed and authenticated
set -e

echo "Cleaning up incorrect 1.0.x releases and tags..."
echo ""

# Check if gh CLI is available
if ! command -v gh &> /dev/null; then
  echo "ERROR: gh CLI is not installed."
  echo "Please install it from: https://cli.github.com/"
  exit 1
fi

# Check if authenticated
if ! gh auth status &>/dev/null; then
  echo "ERROR: Not authenticated with GitHub CLI."
  echo "Please run: gh auth login"
  exit 1
fi

# List of tags to remove
TAGS_TO_REMOVE="1.0.0 1.0.1 1.0.2 1.0.3 1.0.4 1.0.5 1.0.6 v1.0.0"

# Delete GitHub releases (using --cleanup-tag to also remove tags)
for tag in $TAGS_TO_REMOVE; do
  echo "Processing: $tag"
  if gh release view "$tag" &>/dev/null; then
    # --cleanup-tag removes both the release and the tag
    gh release delete "$tag" --yes --cleanup-tag
    echo "  ✓ Deleted release and tag: $tag"
  else
    echo "  - Release $tag not found (already deleted)"
  fi
done

# Delete local tags if they still exist
echo ""
echo "Cleaning up local tags..."
for tag in $TAGS_TO_REMOVE; do
  if git rev-parse "$tag" >/dev/null 2>&1; then
    git tag -d "$tag"
    echo "  ✓ Deleted local tag: $tag"
  fi
done

echo ""
echo "✓ Cleanup complete! Only 0.99.x releases should remain."
echo "Verify at: https://github.com/stanstrup/tidyXCMS/releases"
