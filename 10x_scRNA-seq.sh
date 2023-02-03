#!/bin/bash
# hushu
# 10X genomics scRNA-seq pipeline cellranger

# 10x Software Downloads----Cell Ranger - 7.0.1 (August 18, 2022)
wget -O cellranger-7.0.1.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-7.0.1.tar.gz?Expires=1668014582&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cHM6Ly9jZi4xMHhnZW5vbWljcy5jb20vcmVsZWFzZXMvY2VsbC1leHAvY2VsbHJhbmdlci03LjAuMS50YXIuZ3oiLCJDb25kaXRpb24iOnsiRGF0ZUxlc3NUaGFuIjp7IkFXUzpFcG9jaFRpbWUiOjE2NjgwMTQ1ODJ9fX1dfQ__&Signature=NKQGKy9xMIAerBV0YbA3QLjdQV82Lwlq3QvsK5Koq2QjKFlb2~dsVmDxgi5iUZvwhQXa5PY-glDZgBJHlNGhrrjUdVet5xb7TmgIHRZVbrRyh~VmtG4p1nueaqUXEG5Ht9kOnqIDTfD3bEc9ea-x4Qwvm05qyp3WrV0b3nWdAfUOU2Az2hGcXhixt70C~WGdsqhZ7s57S-aclfMh7fKqiLmqB~vhcbfC7hOkULxqV2iNQGO2orNR7hbGaDHefu1f-YUh5i8kED~I2Fes3iQmnffZAC70tf0NxhSqUa4PBHo~sniy4p4ci5P3PKXcGCHGNUd~KOTlkgK-V6K3ht4gfA__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA"
tar -xzvf cellranger-7.0.1.tar.gz

# download Human reference (GRCh38) dataset required for Cell Ranger.
wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz
tar -xzvf refdata-gex-GRCh38-2020-A.tar.gz
ln -s /home/hushu/hushu/cellranger_scRNA-seq/cellranger-7.0.1/cellrange /usr/bin 
# export PATH=/opt/cellranger-7.0.1:$PATH

# Site check script
cellranger sitecheck > sitecheck.txt
# user can upload their own email to get the audit information to ensure 
# that Cellranger will run smoothly after you generate your own Chromium data
cellranger upload s2051418@ed.ac.uk sitecheck.txt 

cellranger count \
--id=TeTE-03BXY-S-T-scR-5-results \ # a unique run id and output folder name
--transcriptome=/home/hushu/hushu/cellranger_scRNA-seq/refdata-gex-GRCh38-2020-A \
# Path of folder containing 10x-compatible transcriptome reference
--fastqs=/home/hushu/hushu/cellranger_scRNA-seq/scranseqdata/TeTE-03BXY-S-T-scR-5 \
# Path to input FASTQ data
--sample=TeTE-03BXY-S-T-scR-5 \
# Prefix of the filenames of FASTQs to select
--localmem=200 \
# Set max GB the pipeline may request at one time. Only applies to local jobs
--localcores=12
# Set max cores the pipeline may request at one time. Only applies to local jobs
