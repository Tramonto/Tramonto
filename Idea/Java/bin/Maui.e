#!/bin/sh

# Run this script with zero or one arguments: 
# -- if no arguments, MauiDesktop will use the MauiConfig.xml
#    file in etc for Maui configuration
#
# -- if one argument, the argument should be the name of a 
#    Maui configuration file to use in place of etc/MauiConfig.xml
#
# Start up the Maui graphical user interface by running this script


java -DIDEA_HOME=$IDEA_HOME -DHOME=$HOME -Xmx256m -classpath  .:$IDEA_HOME/Idea/Java/classes/idea.jar:$IDEA_HOME/Idea/Java/classes/xerces.jar Maui.Interface.MauiDesktop $1 $2

