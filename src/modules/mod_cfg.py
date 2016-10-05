'''
Created on 23/01/2015

@author: adelpozo
'''


#!/usr/bin/python

import sys, re, os, string, urllib, time, math, random

import  ConfigParser

#######################################################################

def read_cfg_file(cfg_filename):
    
    fi = open(cfg_filename,'r')
    
    config = ConfigParser.ConfigParser()
    config.readfp(fi)
    
    hash_cfg = {}
        
    for field in config.options('INPUT'):
        hash_cfg[field] = config.get('INPUT',field)
    
    for field in config.options('OUTPUT'):
        hash_cfg[field] = config.get('OUTPUT',field)
        
    for field in config.options('SOFTWARE'):
        hash_cfg[field] = config.get('SOFTWARE',field)
    
    for field in config.options('REFERENCE'):
        hash_cfg[field] = config.get('REFERENCE',field)
     
    if config.has_section('CONTROLS'):
        hash_cfg['controls'] = config.items('CONTROLS')
        
    fi.close()
    
    return hash_cfg
