#!/bin/bash
hostname
date
cd /nfs/dust/cms/user/mkomm/HNL/histo/limits
export PATH=/nfs/dust/cms/user/mkomm/HNL/histo/env/bin:$PATH
source activate hnl
which python
which conda
python -u make_hists.py "$@" || return 1
date

