#!/bin/bash
# prepare-news.sh
set -e

# Get parameters from semantic-release
SEMANTIC_VERSION=${SEMANTIC_RELEASE_NEXT_RELEASE_VERSION:-$1}
RELEASE_NOTES=$2

if [ -z "$SEMANTIC_VERSION" ]; then
  echo "ERROR: Next release version not set."
  exit 1
fi

echo "Semantic-release version: $SEMANTIC_VERSION"

# Force version to stay in 0.99.x range for Bioconductor development
# Extract the patch version and increment within 0.99.x
if [[ "$SEMANTIC_VERSION" =~ ^0\.99\.([0-9]+)$ ]]; then
  # Already in 0.99.x range, use as-is
  BIOC_VERSION="$SEMANTIC_VERSION"
elif [[ "$SEMANTIC_VERSION" =~ ^[1-9].*$ ]]; then
  # Version >= 1.0.0, map back to 0.99.x
  # Get current version from DESCRIPTION
  CURRENT_VERSION=$(grep "^Version:" DESCRIPTION | sed 's/Version: //')
  if [[ "$CURRENT_VERSION" =~ ^0\.99\.([0-9]+)$ ]]; then
    PATCH="${BASH_REMATCH[1]}"
    BIOC_VERSION="0.99.$((PATCH + 1))"
  else
    # If current version is not 0.99.x, start at 0.99.1
    BIOC_VERSION="0.99.1"
  fi
else
  # Default case
  BIOC_VERSION="0.99.1"
fi

echo "Preparing NEWS.md and DESCRIPTION for version $BIOC_VERSION..."

# Export for potential use by other scripts
echo "$BIOC_VERSION" > .bioc_version

# Format NEWS.md for R/pkgdown
sed -i 's/^# \[\([0-9]\+\.[0-9]\+\.[0-9]\+\)\].*/## Changes in v\1/' NEWS.md
sed -i 's/^## \[\([0-9]\+\.[0-9]\+\.[0-9]\+\)\].*/## Changes in v\1/' NEWS.md
sed -i 's/^# \([0-9]\+\.[0-9]\+\.[0-9]\+\).*/## Changes in v\1/' NEWS.md
sed -i 's/^## \([0-9]\+\.[0-9]\+\.[0-9]\+\).*/## Changes in v\1/' NEWS.md
sed -i '/^# tidyXCMS/d' NEWS.md
sed -i 's/(\([0-9]\{4\}-[0-9]\{2\}-[0-9]\{2\}\))$/ (\1)/' NEWS.md
sed -i 's/### /### /' NEWS.md
sed -i 's/\[compare\/v[0-9].*//' NEWS.md

# Replace ONLY the first (most recent) version number in NEWS.md with BIOC_VERSION
# This ensures the latest release uses 0.99.x format without touching historical entries
sed -i "0,/## Changes in v[0-9]\+\.[0-9]\+\.[0-9]\+/{s/## Changes in v[0-9]\+\.[0-9]\+\.[0-9]\+/## Changes in v$BIOC_VERSION/}" NEWS.md

# Add commit SHA to the version header for traceability
COMMIT_SHA=$(git rev-parse --short HEAD)
sed -i "0,/## Changes in v$BIOC_VERSION/{s/## Changes in v$BIOC_VERSION/## Changes in v$BIOC_VERSION (commit: $COMMIT_SHA)/}" NEWS.md

# Update DESCRIPTION version
sed -i "s/^Version: .*/Version: $BIOC_VERSION/" DESCRIPTION

# Create git commit with Bioconductor version
echo "Creating git commit for version $BIOC_VERSION..."

# Configure git if not already configured
git config user.name "github-actions[bot]" || true
git config user.email "github-actions[bot]@users.noreply.github.com" || true

# Add modified files
git add NEWS.md DESCRIPTION .bioc_version

# Create commit with Bioconductor version (not semantic version)
git commit -m "chore(release): $BIOC_VERSION [skip ci]

$RELEASE_NOTES"

# Push the commit using GITHUB_TOKEN
# The git remote is already configured by semantic-release checkout
git push

echo "Committed and pushed version $BIOC_VERSION"
