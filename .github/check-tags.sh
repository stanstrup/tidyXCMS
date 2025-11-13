#!/bin/bash
# check-tags.sh - Check what tags exist on remote
set -e

echo "Fetching tags from remote..."
git fetch --tags --prune origin

echo ""
echo "=== Tags on remote (1.x.x, 2.x.x, semantic-release-*) ==="
git ls-remote --tags origin | grep -E "(refs/tags/(v?[12]\.|semantic-release-))" || echo "No incorrect tags found!"

echo ""
echo "=== Local tags (1.x.x, 2.x.x, semantic-release-*) ==="
git tag | grep -E "^(v?[12]\.|semantic-release-)" || echo "No incorrect tags found!"

echo ""
echo "=== All remote tags ==="
git ls-remote --tags origin | cut -f2 | sed 's|refs/tags/||' | sort -V
