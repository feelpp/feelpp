#!/bin/bash

# Determine the workspace directory
WORKSPACE_DIR="${1:-/workspaces/feelpp}"

# Define the file path
FILE_PATH="/home/feelpp/.feelppconfig"

# Create the JSON content with the dynamic workspace directory
cat > "$FILE_PATH" <<EOF
{
 "append_date": false,
 "append_np": true,
 "feelppdb": "feelppdb",
 "exprs": "exprs",
 "geos": "geo",
 "location": "global",
 "logs": "logs",
 "owner": {
  "email": "",
  "name": "feelpp"
 },
 "global_root":"$WORKSPACE_DIR/feelppdb"
}
EOF

echo "File created at $FILE_PATH, global_root=$WORKSPACE_DIR/feelppdb"
