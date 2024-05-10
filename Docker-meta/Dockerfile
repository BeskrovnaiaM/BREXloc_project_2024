FROM python:3.12

WORKDIR /home

COPY requirements.txt requirements.txt

RUN pip install -r requirements.txt

COPY brextractor.py brextractor.py

RUN mkdir brextraction

ENTRYPOINT bash

