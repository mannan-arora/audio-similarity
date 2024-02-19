#!/bin/bash

# Run make
make

# Check if make was successful
if [ $? -eq 0 ]; then
    # Run the executable
    ./similarity
    
    # Delete object & depenency files
    rm *.o
    rm *.d
fi