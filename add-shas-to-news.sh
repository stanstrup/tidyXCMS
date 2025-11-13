#!/bin/bash
# add-shas-to-news.sh
# This script adds commit SHAs to all version headers in NEWS.md
# Run this once to add traceability to all historical releases

set -e

echo "Adding commit SHAs to all version headers in NEWS.md..."

# Create a temporary file
cp NEWS.md NEWS.md.backup

# Get all release commits and their versions
git log --grep="chore(release):" --all --format="%H" | while read commit; do
  version=$(git show $commit:DESCRIPTION 2>/dev/null | grep "^Version:" | sed 's/Version: //')
  if [ -n "$version" ]; then
    short_sha=$(git rev-parse --short $commit)
    echo "  v$version -> commit $short_sha"

    # Replace the version header if it doesn't already have a commit SHA
    sed -i "s/^## Changes in v$version\$/## Changes in v$version (commit: $short_sha)/" NEWS.md
  fi
done

echo ""
echo "✓ All version headers updated with commit SHAs"
echo ""
echo "Backup saved to NEWS.md.backup"
echo "Review changes with: git diff NEWS.md"
