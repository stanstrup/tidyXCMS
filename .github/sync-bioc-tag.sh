#!/bin/bash
# sync-bioc-tag.sh - Create 0.99.x tag and GitHub release
set -e

# Get parameters from semantic-release
SEMANTIC_VERSION=$1
RELEASE_NOTES=$2

echo "Semantic-release version: $SEMANTIC_VERSION"

# Read the Bioconductor version that was set by prepare-news.sh
if [ ! -f ".bioc_version" ]; then
  echo "ERROR: .bioc_version file not found"
  exit 1
fi

BIOC_VERSION=$(cat .bioc_version)
echo "Creating Bioconductor release: $BIOC_VERSION"

# Create lightweight tag for the Bioconductor version
git tag -f "$BIOC_VERSION"

# Push the tag to remote
git push origin "$BIOC_VERSION" --force

echo "Created and pushed tag $BIOC_VERSION"

# Extract release notes from NEWS.md (first section after the version header)
# This ensures we use the properly formatted NEWS.md content
NEWS_CONTENT=$(awk '/## Changes in v'"$BIOC_VERSION"'/,/## Changes in v[0-9]/ {
  if (/## Changes in v[0-9]/ && !/## Changes in v'"$BIOC_VERSION"'/) exit;
  if (!/## Changes in v'"$BIOC_VERSION"'/) print
}' NEWS.md)

# Create GitHub release using gh CLI
if [ -n "$GITHUB_TOKEN" ]; then
  echo "Creating GitHub release for $BIOC_VERSION..."

  # Check if release already exists
  if gh release view "$BIOC_VERSION" &>/dev/null; then
    echo "Release $BIOC_VERSION already exists, deleting it first..."
    gh release delete "$BIOC_VERSION" --yes
  fi

  # Create the release with NEWS.md content
  gh release create "$BIOC_VERSION" \
    --title "tidyXCMS v$BIOC_VERSION" \
    --notes "$NEWS_CONTENT" \
    --latest

  echo "GitHub release created successfully for $BIOC_VERSION"
else
  echo "WARNING: GITHUB_TOKEN not set, skipping GitHub release creation"
fi
