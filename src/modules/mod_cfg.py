#!/usr/bin/
"""
Created on 23/01/2015

@author: adelpozo
"""


import sys
import re
import os
import string
import urllib
import time
import math
import random

import ConfigParser

#######################################################################


def read_cfg_file(cfg_filename):
    """

    :param cfg_filename:
    :return:
    """
    fi = open(cfg_filename, 'r')
    
    config = ConfigParser.ConfigParser()
    config.readfp(fi)
    
    hash_cfg = {}
    for field in ['INPUT', 'OUTPUT', 'SOFTWARE', 'REFERENCE']:
        for field_value in config.options(field):
            hash_cfg[field_value] = config.get(field, field_value)

    if config.has_section('CONTROLS'):
        hash_cfg['controls'] = config.items('CONTROLS')

    fi.close()
    
    return hash_cfg
