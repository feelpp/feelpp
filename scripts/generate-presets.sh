#!/bin/bash
# Generate CMakePresets.json from CMakePresets.jsonnet
# Usage: ./scripts/generate-presets.sh

set -e

cd "$(dirname "$0")/.."

# Check if jsonnet is installed
if ! command -v jsonnet &> /dev/null; then
    echo "Error: jsonnet is not installed"
    echo "Install with:"
    echo "  Ubuntu/Debian: sudo apt install jsonnet"
    echo "  macOS: brew install jsonnet"
    echo "  Or via Go: go install github.com/google/go-jsonnet/cmd/jsonnet@latest"
    exit 1
fi

echo "Generating CMakePresets.json from CMakePresets.jsonnet..."
jsonnet CMakePresets.jsonnet | python3 -m json.tool --indent 3 --no-ensure-ascii > CMakePresets.json

echo "✓ Generated CMakePresets.json successfully"
echo "  Original presets: $(grep -c '"name":' CMakePresets.json.backup 2>/dev/null || echo 'N/A')"
echo "  Generated presets: $(grep -c '"name":' CMakePresets.json)"
