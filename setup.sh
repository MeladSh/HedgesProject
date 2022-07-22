#!/bin/bash

install_requirements(){
  echo "--------------------Installing Python3.8--------------------"
  sleep 1
  sudo apt-get install python3.8-dev
  sleep 1
  echo "--------------------Installing numpy-------------------"
  sleep 1
  sudo apt-get install python3.8-numpy
}

run(){
  # install requirements -> Python3.8 / numpy
  install_requirements
}

run
