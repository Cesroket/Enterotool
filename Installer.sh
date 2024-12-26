#!/bin/bash

echo "Installing dependencies"
#the idea is to make Enterotool as a variable in .bashrc

## this is the function that will be added to .bashrc
function Enterotool() {

if [[ "$1" == "--mkdb" ]]; then
    echo "making database"
    python ~/Documents/Enterobaser.sh

elif [[ "$1" == "--h" ]]; then
    cat /home/cesar/Documents/Enterotool-manual.txt
else
    bash ~/Documents/Enterotool.sh "$@"

fi
}
