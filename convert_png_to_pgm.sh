#!/bin/bash

# Check if a folder name is provided as a parameter
if [ -z "$1" ]; then
    echo "Usage: $0 <folder_name>"
    echo "Example: $0 video or $0 my_video"
    exit 1
fi

# Set the folder name based on the parameter
VIDEO_DIR="$1"

# Check if the specified folder exists
if [ ! -d "$VIDEO_DIR" ]; then
    echo "Error: Folder '$VIDEO_DIR' not found"
    exit 1
fi

# Loop through all png files in the specified folder
for png_file in "$VIDEO_DIR"/*.png; do
    # Get the filename without extension
    filename=$(basename "$png_file" .png)
    
    # Convert png to PGM and save in the output folder
    convert "$png_file"  "${VIDEO_DIR}/${filename}.pgm"
    
    echo "Converted $png_file to ${VIDEO_DIR}/${filename}.pgm"
done

echo "Conversion complete. PGM files are in the '${VIDEO_DIR}' folder."