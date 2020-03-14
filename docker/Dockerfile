FROM dealii/dealii:latest

MAINTAINER luca.heltai@gmail.com

USER root

RUN apt-get update && apt-get -yq install \
    vim \
    emacs-nox \
    gdb

USER dealii
WORKDIR /home/dealii
