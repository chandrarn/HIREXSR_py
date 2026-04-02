
# HIREXSR_py
A python port of THACO IDL analysis scripts for the Alcator C-Mod HIREX Doppler spectrometry data.
The original files can be found in: /usr/local/cmod/idl/HIREXSR/

## Installation
The project is set up to run with uv: use uv sync to drop into the correct environment, and then uv run for the codes of your choice

## Usage
See the jupyter notebook for run examples on pulling and plotting the line integrated or inverted data, or checking a list of shots to see if good data exists.

## Features
The code features a number of switches to account for data being missing on some MDS nodes, or nodes potentially being renamed.
Users may need to manually check a number of lines to find good data

## Configuration
Project is set up to run within a uv envionrment. If the user does not wish to do this, the only potentially nonstanard package is mdsthin, the python port of the original MDS+ code for data access.
