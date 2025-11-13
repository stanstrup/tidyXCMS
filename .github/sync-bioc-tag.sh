#!/bin/bash
# sync-bioc-tag.sh - Create and sync 0.99.x tag with GitHub release
set -e

# Read the Bioconductor version that was set by prepare-news.sh
if [ -f ".bioc_version" ]; then
  BIOC_VERSION=$(cat .bioc_version)
  echo "Syncing Bioconductor version tag: $BIOC_VERSION"

  # Create lightweight tag for the Bioconductor version
  git tag -f "$BIOC_VERSION"

  # Push the tag to remote (requires authentication via GITHUB_TOKEN)
  git push origin "$BIOC_VERSION" --force

  echo "Created and pushed tag $BIOC_VERSION"
else
  echo "WARNING: .bioc_version file not found, skipping tag sync"
fi
