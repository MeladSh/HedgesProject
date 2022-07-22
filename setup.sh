#!/bin/bash

install_requirements(){
  echo "--------------------Installing Python3.8--------------------"
  sleep 1
  sudo apt-get install python3.8-dev
  sleep 1
  echo "--------------------Installing numpy1.17.4-------------------"
  sleep 1
  sudo apt-get install python3.8-numpy
}

run(){
  # install requirements -> Python3.8 / numpy1.17.4
  install_requirements
}

run
