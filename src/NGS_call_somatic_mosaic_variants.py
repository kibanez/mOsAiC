#!/usr/bin/python

import sys, re, shlex , os, string, urllib, time, math, random, subprocess, shutil
from os import path as osp
import optparse
#from wx.lib.agw.thumbnailctrl import file_broken
from HTSeq import SAM_Reader
import ConfigParser
from subprocess import Popen, PIPE
from itertools import izip_longest, chain
import vcf
from vcf.parser import _Filter,_Record,_Info,_Format,_SV
from vcf.utils import walk_together
from operator import itemgetter
import pandas as pd
from tempfile import mkstemp
from shutil import move
from os import remove, close
import warnings
import numpy as np
import logging
from modules import mod_cfg
from modules import mod_variant
#import xlsxwriter

modulesRelDirs = ["modules/"]

localModulesBase = osp.dirname(osp.realpath(__file__))

for moduleRelDir in modulesRelDirs:
        sys.path.insert(0,osp.join(localModulesBase, moduleRelDir))



principalsRelDirs = ["../"]
for pRelDir in principalsRelDirs:
        sys.path.insert(0,osp.join(localModulesBase, pRelDir))


class OptionParser(optparse.OptionParser):

    def check_required(self, opt):

        option = self.get_option(opt)

        atrib = getattr(self.values, option.dest)
        
        if atrib is None:
            self.error("%s option not supplied" % option)

########################################################################################################################

def change_sample_name(vcf_file,new_label='overlapping'):
    
    vcf_reader = vcf.Reader(filename=vcf_file)
    
    l_records = []
    
    for record in vcf_reader:
        record.samples = [record.samples[0]]
        l_records.append(record)
        
    vcf_reader.samples = [new_label]
    
    vcf_tmp = os.path.splitext(vcf_file)[0]+'.tmp.vcf'
    
    vcf_writer = vcf.Writer(open(vcf_tmp , 'w'), vcf_reader)
    
    for record in l_records:
        vcf_writer.write_record(record)
    
    return vcf_tmp

########################################################################################################################

def combine_vcf_v2(l_vcf_files,l_rod_priority,mode,ouput_path=None,**kwargs):
    
    """
        There are two modes:
        
        1. merge: this is when you need to join hetereogeneous calls into the same register with the condition that they overlap into the same genomic position. 
                    It could be merged calls from different tools belonging to the same sample
        2. combine: When you need to join the same variants (shared ref and alt alleles) into one register. Its recommended in familiar studies. In this mode, only different samples are allowed 
    """
        
    if mode == "merge":
        get_key = lambda r: (r.CHROM, r.POS)
    elif mode == "combine":
        get_key = lambda r: (r.CHROM, r.POS, r.REF, r.ALT)
    else:
        raise InputError('mod_variant.combine_vcf_v2: The mode of combining vcf is not correct: %s' % (mode))
    
    ### specific configurations
    label_mode = "no_append" ### possibilities: append/no_append
    if kwargs.has_key('label_mode'):
        label_mode = kwargs['label_mode']
    
    info_mode = "overwrite" ### possibilities: overwrite/append
    if kwargs.has_key('info_mode'):  
        info_mode = kwargs['info_mode']
        
    remove_non_variants = False ### possibilities: True/False
    if kwargs.has_key('remove_non_variants'):
        remove_non_variants = kwargs['remove_non_variants'] 
        
    ### 1. generate the name of the combined file    
    vcf_file_combined = ""
        
    root_name = __get_root_name(l_vcf_files)
    
    if root_name == "":
        root_name = "v"
    
    if ouput_path == None:
        vcf_file_combined = root_name+'.combined.vcf'
    else:
        vcf_file_combined = os.path.join(ouput_path,root_name+'.combined.vcf')
    
    
    ### 2. Check out the arragement of the vcf files in function of the rod list 
    if l_rod_priority == [] or l_rod_priority == None:
        l_rod_priority = l_vcf_files
        label_mode = "no_append"
        
    if len(l_rod_priority) != len(l_vcf_files):
        raise InputError('mod_variant.combine_vcf_v2: The number of ids in the list of rod priority do not agree with the number of vcf files')
    
    num_files = len(l_vcf_files)
    
    header = ""
    
    
    ### 3. Fill the list of vcf readers. In addition, it configure the INFO field in the header
    ### Also, it assigns the index of the file and the samples
    ### It writes the header into a string 
    l_readers = []
    
    hash_vcf_2_readers = {}
    hash_sample_2_file = {}
    hash_sample_2_index= {}
    hash_sample_2_rod  = {}
    
    i_sample = 0
    
    dict_filters = {}
    
    for vcf_file in l_vcf_files:
        r = vcf.Reader(filename=vcf_file)
        dict_filters.update(r.filters)
    
    for (i,(vcf_file,rod)) in enumerate(zip(l_vcf_files,l_rod_priority)):
        
        r = vcf.Reader(filename=vcf_file)
        
        r.infos['NC']  = _Info('NC',1,'Integer',"Number of Calls")
        r.infos['NV']  = _Info('NV',1,'Integer',"Number of Variants")
        r.infos['set'] = _Info('set',1,'String',"Source of the vcf file")
                        
        r.filters = dict_filters
        
        l_readers.append(r)
        
        if header == "":
            header = __get_header(r)
        if header.find('##contig') == -1:
            try:
                hash_contigs = r.contigs
                for cont in hash_contigs:
                    header += "##contig=<ID=%s,length=%d>\n" % (hash_contigs[cont].id,hash_contigs[cont].length)  
            except AttributeError:
                pass
                        
        hash_vcf_2_readers[os.path.basename(vcf_file)] = i
        
        l_s = r.samples # samples of each vcf
        
        for sample in l_s:
            
            if mode == "combine":
                if hash_sample_2_file.has_key(sample):
                    raise RuntimeError("combine_vcf:combine_vcf_v2: There are at least two vcf files with the same sample label\n")
            
            if mode == "merge": # the sample could be in two o more vcf files
                hash_sample_2_file.setdefault(sample,[])
                hash_sample_2_file[sample].append(os.path.basename(vcf_file))
                
                hash_sample_2_rod.setdefault(sample,[])
                hash_sample_2_rod[sample].append(rod)
            else:
                hash_sample_2_file[sample] = os.path.basename(vcf_file)
                hash_sample_2_rod[sample]  = rod
            
            if not hash_sample_2_index.has_key(sample): # in case of duplicated samples, it assigns the lowest index
                hash_sample_2_index[sample] = i_sample
                i_sample += 1
            
    l_total_samples =  map(itemgetter(1),sorted(map(lambda (i,s): (i,s), zip(hash_sample_2_index.values(),hash_sample_2_index.keys()))))
    
    num_samples = len(l_total_samples)
    
    if mode != "combine":

        fo = open(vcf_file_combined,'w')
        fo.write(header)
        fo.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n" % ('\t'.join(l_total_samples)))
    
    ### 4. It iterates over each reader (vcf_file) and then, for the samples of the reader
    elif mode == "combine":
        hash_walker = {}
    
    for l_record in vcf.utils.walk_together(*l_readers,vcf_record_sort_key=get_key):
       
        sum_QUAL = 0
        num_r    = 0
        num_s    = 0
        filter_string   = ""
        set_str  = ""
        is_PASS  = False 
        
        hash_format = {}
        
        ### It counts the number of calls. A register with None means no call ./.
        num_calls = 0
        ### It counts the number of variants. A register with 0/0 is not counted
        num_variants = 0
        ### It counts the number of vcf with calls
        num_records = 0
        for r in l_record:
            if (type(r) == vcf.model._Record):
                num_records += 1
        
        l_format  = ['GT','AD','DP','GQ','PL'] # these are the fields of the genotype taken into account 
        # by default, the no-call is represented so
        hash_none = {'GT':'./.'} 
        l_gt      = num_samples*[['./.','0,0','0','.','.']]
        
        chrm     = ""
        pos      = 0
        ref      = 0
        l_alt    = []
        l_id     = []
        l_filter = []
        l_set = map(lambda k: None, range(num_samples))
        
        hash_filter = {}
        hash_info_all = {}
        
        for r_i in l_record:
            
            if type(r_i) != vcf.model._Record:
                continue
            
            chrm  = r_i.CHROM
            pos   = r_i.POS
            ref   = r_i.REF
            l_alt_r = map(lambda sb: sb.__repr__() , r_i.ALT)
            #l_s     = map(lambda s: s.sample, r_i.samples)
            
            l_s = []
            for c in r_i.samples:
                if c.gt_bases != None:
                    l_s.append(c.sample)
            
            if l_alt_r == ['None']: # When the site has being for all the samples genotyped and is HOM_REF
                l_alt_r   = ['.']
                 
            l_alt.extend(l_alt_r)
            
            ### Build the ALT fields
            l_alt = list(set(l_alt))
            alt_str = ','.join(l_alt)
            if l_alt != ['.']:
                if set(l_alt).intersection(set(['.']))!= set([]):
                #if filter(lambda a: a=='.', l_alt) != []:
                    alt_str = alt_str.replace('.', '')
                    alt_str = alt_str.replace(',', '') 
            
                            
            hash_info = dict(r_i.INFO)
            
            is_HOM_REF = False
            is_MULTIALLELIC = False
            
            for sample in l_s:
                if r_i.genotype(sample)['GT'] == "0/0":
                    is_HOM_REF = True
                    num_calls += 1
                elif r_i.genotype(sample)['GT'] != "0/0" and r_i.genotype(sample)['GT'] != None:
                    num_calls    += 1
                    num_variants += 1
                     
            if len(l_alt_r) > 1:
                is_MULTIALLELIC = True
            
            if r_i.ID != None:
                l_id.extend(r_i.ID.split(';'))
            else:
                if not is_HOM_REF:
                    l_id.append('.')
                elif is_MULTIALLELIC:
                    l_id.append('.')
                else:
                    l_id = ['.']
            
            if r_i.QUAL != None:
                sum_QUAL += r_i.QUAL
            
            if r_i.QUAL != None:
                num_r += 1
            
            for sample in l_s: # Iterating across the samples of a vcf file
            
                i_sample = hash_sample_2_index[sample] # this index is unique 
                
                if mode == "combine":
                    
                    if label_mode == "no_append":
                        # Instead of printing the absolute path, only the DNA name
                        #rod   = hash_sample_2_rod[sample]
                        str_aux = sample.split('-')
                        rod   = str_aux[len(str_aux)-1]
                        
                    else: # mode append
                        rod   = hash_sample_2_rod[sample] + '_%s' % (sample)
                        
                    f_v   = hash_sample_2_file[sample] # this only to track the files that the sample is coming from
                    i_vcf = hash_vcf_2_readers[f_v] # this only to track the position in the rod priority list of the file/s that the sample is coming from
                else: ### revisar!!!!!
                    if label_mode == "no_append":
                        rod   = ','.join(hash_sample_2_rod[sample])
                    else: # mode append
                        rod   = ','.join(map(lambda r: r+ "_%s" % (sample), hash_sample_2_rod[sample]))
                    l_f_v = hash_sample_2_file[sample]
                    l_i_vcf = map(lambda f_v: hash_vcf_2_readers[f_v], l_f_v)
                            
                hash_gt = {} # hash of the genotipe information
                
                if r_i.FILTER != [] and r_i.FILTER != None:
                    l_filter.extend(r_i.FILTER) 
                    l_set[i_sample] = "filterIn"+rod
                    hash_filter[sample] = True
                else:
                    l_set[i_sample] = rod
                    is_PASS = True
                    
                if info_mode == "append":
                   hash_info_all.update(hash_info) 
            
                try:
                    hash_gt['GT'] = r_i.genotype(sample)['GT']
                except:
                    pass
                
                try:
                    hash_gt['DP'] = str(r_i.genotype(sample)['DP'])
                except:
                    if is_HOM_REF: 
                        dp = hash_info.get('DP',[0])
                        if type(dp) == int:
                            hash_gt['DP'] = str(dp)
                        else:
                            hash_gt['DP'] = str(dp[0])
                    pass
            
                try:
                    AD = r_i.genotype(sample)['AD']
                    if type(AD) == int:
                        AD = (int(hash_gt['DP'])-AD,AD)
                    hash_gt['AD'] = ','.join(map(str,AD))
                except:
                    pass
            
                try:
                    hash_gt['GQ'] = str(r_i.genotype(sample)['GQ'])
                except:
                    if is_HOM_REF: 
                        hash_gt['GT'] = "0/0"
                    pass
            
                
                try:
                    hash_gt['PL'] = ','.join(map(lambda p: str(p).replace('None','.'), r_i.genotype(sample)['PL']))
                except:
                    pass
                
                l_ = []
                for ft in l_format:
                    ft_gt = hash_gt.get(ft,'.')
                    if ft == "GT" and str(ft_gt) == 'None':
                         ft_gt = "./."
                    elif ft == "AD" and (ft_gt == '.' or str(ft_gt) == 'None'):
                        ft_gt = "0,0"
                    elif ft == "DP" and (ft_gt == '.' or str(ft_gt) == 'None'):
                        ft_gt = "0"
                    elif ft == "PL" and (ft_gt == '.' or str(ft_gt) == 'None'):
                        ft_gt = "."
                    elif ft == "GQ" and (ft_gt == '.' or str(ft_gt) == 'None'):
                        ft_gt = '.'
                    l_.append(ft_gt)
            
                l_gt[i_sample] = l_
       
       
        # We check here whether the walk together generator has set in different lines equal key elements (we have seen that it does...!!)
        if mode == "combine" and hash_walker.has_key((chrm,pos,ref,alt_str)):
            key_repeated = True
        else:
            key_repeated = False
        
        ### Mean quality over the records with qual value assigned
        QUAL = 0
        if (num_r >= 1):        
            QUAL = float(sum_QUAL)/num_r
        
        ### Generating the strings of the samples which calls overlap.
        ### Intersection -> means that all the samples (from all the vcf files) have variants and share the same context (i.e. merge or combine context)
        l_set   = sorted(list(set(l_set)))
        
        set_str = ""
        if (num_variants == num_samples):
            set_str = "Intersection"
        else:
            for s in l_set:
                if s!=None:
                    set_str += s+'--'
            set_str = set_str[:-2]
        
                
        ### Generating the string of the FILTER field. If any of the records with call is filtered, it is included in the FILTER string
        ### PASS means all the calls do PASS
        num_filtered = len(hash_filter.keys())
        l_filter = list(set(l_filter))
        filter_string  = "PASS"
        if l_filter != []: # in case any call is filtered
            #if num_calls == num_filtered:
            # JUNE'16: OJO!!! if all the samples in the pool are filtered, so, filteredInALL and NC = num_samples. Otherwise, we do not know which samples are filtered
            if num_samples == num_filtered:
                set_str = "FilteredInAll"
                
            filter_string = ""
            for f in l_filter:
                if f!='PASS':
                    filter_string += f+';'  
            filter_string = filter_string[:-1]
        

        # we check whether the sample in set is wthin the set list already saved in the hash_walker            
        # the set field which inclued all the samples having the variant        
        if key_repeated:
        
            l_set_walker = hash_walker[(chrm,pos,ref,alt_str)]['info'].split(';')[0]
            l_set_walker = l_set_walker.split('--')            
            
            num_samples_walker = int(hash_walker[(chrm,pos,ref,alt_str)]['info'].split(';')[1].split("=")[1])
            
            
            # l_set_walker is a list as follows: ['set=*vcf','*vcf',...,filterInvcf*]
            l_set_walker = clean_walker(l_set_walker)     
            
            # sometimes l_set is full of None's because the vcf module cannot found in the list of alleles the allele number
            
            # l_set is the list of the samples with the variant (the first element is always a None)
            if l_set[0] == None:
                l_set = l_set[1:]
            
            if len(l_set) > 0:                
                        
                l_set_unique = str(l_set[0])
                
                # number of samples that were not included before in the same key hash_walker (in order to include them now)
                # this is FUNDAMENTAL since "set=FilteredInAll" is even when NC=1, NC=12, etc. 
                new_num_samples = 0
                
                if len(l_set) > 1:
                    # check whether rod is a list
                    for s in l_set:                
                        if str(s) not in l_set_walker:
                            # we check whether the sample is filter in (previously if it has "filterIn" it does, otherwise no
                            if "filterIn" not in str(s):
                                filter_s = "PASS"
                            else:
                                filter_s = "filtered"
                                
                            if filter_s != 'PASS':
                                l_set_walker.append(s)
                            else:
                                l_set_walker = [s] + l_set_walker
                            
                            new_num_samples = new_num_samples + 1
                else:
                    
                    new_num_samples = 1
                    
                    if l_set_walker == ["FilteredInAll"]:

                        num_calls = num_samples
                    
                    else:
                                 
                        if l_set_unique not in l_set_walker:
                            if "filterIn" not in l_set_unique:
                                filter_s = "PASS"
                            else:
                                filter_s = "filtered"
                           
                            if filter_s != 'PASS' and l_set_walker == ["FilteredInAll"]:
                                l_set_walker = ["FilteredInAll"]
                            elif filter_s != 'PASS' and l_set_walker != ["FilteredInAll"]:
                                l_set_walker.append(l_set_unique)    
                            else:
                                l_set_walker = [l_set_unique] + l_set_walker
            
                set_str = '--'.join(l_set_walker)
                
                if set_str != "FilteredInAll":
                    num_calls = len(l_set_walker)
        
            else:
                
                continue
            
        ## Build the INFO field
        info    = "set=%s;NC=%d;NV=%d" % (set_str,num_calls,num_variants)
        if info_mode == "append":
            str_to_append = ""
            try:
                l_keys = list(set(hash_info_all.keys()).difference(set(['set','NC','NV'])))
                for k in l_keys:
                    value = hash_info_all[k]
                    if type(value) == list:
                        value = ','.join(map(str,value))
                    str_to_append += ";" + "%s=%s" % (k,value)
            except:
                pass
            
            info = info+str_to_append
            if info[-1] == ";":
                info = info[:-1]
                
        ## Build the FORMAT field
        if remove_non_variants:
            if filter(lambda gt: gt == '0/0' or gt == 'None', l_gt) != []:
                continue
        
        gt_all  = map(lambda i: ':'.join(i) , l_gt)
        ## Build the ID field
        l_id = filter(lambda id: id!='.', list(set(l_id)))
        if l_id == []:
            l_id = ['.']
        id = ','.join(list(set(l_id)))
        
#         ### Build the ALT fields
#         l_alt = list(set(l_alt))
#         alt_str = ','.join(l_alt)
#         if l_alt != ['.']:
#             if set(l_alt).intersection(set(['.']))!= set([]):
#             #if filter(lambda a: a=='.', l_alt) != []:
#                 alt_str = alt_str.replace('.', '')
#                 alt_str = alt_str.replace(',', '') 
        
        if mode == "combine":
           
            hash_variant = {}
            
            hash_variant['SNP_id'] = str(id)
            hash_variant['qual'] = str(QUAL)
            hash_variant['filter'] = filter_string
            hash_variant['info'] = info
            hash_variant['format'] = ':'.join(l_format)
            hash_variant['format_info'] = '\t'.join(gt_all)
            
            
            hash_walker[(chrm,pos,ref,alt_str)] = hash_variant
        
        
        else:
            
            #CHROM    POS    ID    REF    ALT    QUAL    FILTER    INFO    FORMAT    40204_Illumina    40204_bowtie    40204_bwa
            fo.write("%s\t%d\t%s\t%s\t%s\t%1.2f\t%s\t%s\t%s\t%s\n" % (chrm,pos,id,ref,alt_str,QUAL,filter_string,info,':'.join(l_format),'\t'.join(gt_all)))    
    
    if mode == "combine":
        
        # check if hash contains new
        fo = open(vcf_file_combined,'w')
        fo.write(header)
        fo.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n" % ('\t'.join(l_total_samples)))
        
        #aqui
        l_cyto_template    = map(lambda i: "chr%i" % (i), range(1,23))+['chrX','chrY','chrM']
        hash_cyto_template = dict(zip(l_cyto_template,range(len(l_cyto_template))))
        hash_cyto          = {}
                
        for chrm in list(set(map(itemgetter(0),hash_walker.keys()))): 
        
            if hash_cyto_template.has_key(chrm):
                hash_cyto[chrm] = hash_cyto_template[chrm] 
            else:
                chrm_2 = chrm.split('_')[0]
                if hash_cyto_template.has_key(chrm_2):
                    hash_cyto[chrm] = hash_cyto_template[chrm_2]
                else:
                    raise AttributeError('mod_variant.combine_v2: There is chromosome entry not found: %s' % (chrm))
                
        chr_cyto_sorted = map(itemgetter(0),sorted(map(lambda c,p: (c,p) , hash_cyto.keys(),hash_cyto.values()),key=itemgetter(1,0)))
        
        l_fields = ['qual','filter','info','format','format_info']
                
        for key in sorted(sorted(hash_walker.keys(),key=itemgetter(1)),key=lambda x: chr_cyto_sorted.index(x[0])):        
            
            l_key = map(str,list(key))
            l_key.insert(2,hash_walker[key].get('SNP_id','.'))
            key_str = '\t'.join(l_key)
            
            fo.write(key_str + '\t' + '\t'.join(map(lambda field: hash_walker[key].get(field,'.'), l_fields))+'\n')
            
        
        fo.close()
    
    else:
    
        fo.close()
            
    return vcf_file_combined

########################################################################################################################
# we create a table with all the variants corresponding ONLY to the tissue variants and NOT to the control
# only if blood is True, will be add to "tejido" when lowdepthmosaicism

def parse_control(vcf_input,hash_tissue,blood=False):
    # hash_table will have, those variants detected only in tissue and not in blood (with PASS filter)
    # only if blood has low_depth ==> we keep it with that filter flag
    
    hash_table = hash_tissue
    
    vcf_reader = vcf.Reader(filename=vcf_input)
    
    for i,r in enumerate(vcf_reader):
        
        try:
            
            if type(r) ==  _Record:
                
                hash_fields = dict(r.INFO)
                #hash_fields.update({'low_mosaic':''})
                
                chr = r.CHROM
                pos = r.POS    
                ref = str(r.REF)
                alt = r.ALT
                
                # ALT
                if alt != [None]:
                    l_alt = []
                    if len(alt) > 1:                        
                        for x in alt:
                            l_alt.append(str(x))
                        alt_str = ','.join(l_alt)
                    else:
                        alt_str = str(alt[0])
                        l_alt = alt_str
                
                alt = alt_str
                

                key = (chr,pos,ref,alt)
                
                if not hash_table.has_key(key):
                    continue  
                         
                # FILTER
                l_filter = []
                
                if r.FILTER != [] and r.FILTER != None:
                    l_filter.extend(r.FILTER) 
                else:
                    l_filter.append("PASS")
                
                
                ### Generating the string of the FILTER field. If any of the records with call is filtered, it is included in the FILTER string
                ### PASS means all the calls do PASS
                filter = str(l_filter[0])
                
                if blood == False:
                    
                    if hash_table.has_key(key):
                        hash_table.pop(key, None)
                                            
                else: # we check the filter. if PASS, we remove it, otherwise we add a new filed: Low_Mosaic
                       
                    if filter == "PASS":
                        if hash_table.has_key(key):
                            hash_table.pop(key, None)
                    else:
                       
                        mavaf = str(hash_fields.get('MAVAF','.'))
                        # change low_mosaic info in hash_table
                        l_info = hash_table[key]['info'].split(";")
                        
                        new_info = []
                        for j in l_info:
                            aux = j.split("=")
                            
                            if (aux[0] == "low_mosaic"):
                                new_data = aux[0] + "=" + mavaf
                            else:
                                new_data = j
                            
                            new_info.append(new_data)
                        
                        info = ';'.join(new_info)
                        
                        # update ONLY the low_mosaic field in info to the mosaic hash
                        hash_table[key]['info'] = info

        except:
            
            raise RuntimeError('mosaic_somatic_variants_calling.parse_control: Some error has occurred in variant line %i:\n%d' % (i))    
    
    
    return hash_table

########################################################################################################################
# we create a table with all the variants corresponding to the tissue variants
def parse_tissue_vcf(vcf_input):
    
    hash_table = {}
    
    vcf_reader = vcf.Reader(filename=vcf_input)
    
    vcf_reader.infos['low_mosaic']  = _Info('low_mosaic',1,'String',"The percentage of the mosaicism detected in the blood vcf")
    
    header = __get_header(vcf_reader)

    for i,r in enumerate(vcf_reader):
        
        try:
            
            if type(r) ==  _Record:
                
                hash_fields = dict(r.INFO)
                hash_fields.update({'low_mosaic':'0'})
                
                chr = r.CHROM
                pos = r.POS
                ref = str(r.REF)
                alt = r.ALT
                qual = r.QUAL
                if qual == None:
                    qual = "."
                
                l_info = []
                for (key,val) in hash_fields.items():
                    if isinstance(hash_fields[key],float):
                        continue
                    if not isinstance(hash_fields[key],str) and not isinstance(hash_fields[key],int):
                        aux = key + "=" + ','.join(map(str,hash_fields[key]))
                    else:
                        aux = key + "=" + str(hash_fields[key])
                        
                    l_info.append(aux)
                #info = ';'.join("{!s}={!r}".format(key,val) for (key,val) in hash_fields.items())
                info = ";".join(l_info)
                format = r.FORMAT
                
                
                
                # ID
                if r.ID == None:
                    id = '.'
                
                
                # ALT
                if alt != [None]:
                    l_alt = []
                    if len(alt) > 1:                        
                        for x in alt:
                            l_alt.append(str(x))
                        alt_str = ','.join(l_alt)
                    else:
                        alt_str = str(alt[0])
                        l_alt = alt_str
                
                alt = alt_str
                
                hash_table[(chr,pos,ref,alt)] = {}
                
                # FILTER
                l_filter = []
                
                if r.FILTER != [] and r.FILTER != None:
                    l_filter.extend(r.FILTER) 
                else:
                    l_filter.append("PASS")
                
                
                ### Generating the string of the FILTER field. If any of the records with call is filtered, it is included in the FILTER string
                ### PASS means all the calls do PASS
                filter = str(l_filter[0])
                
                
                # FORMAT
                sample = r.samples[0].sample
                
                l_format  = ['GT','PL','AD']
                l_ = []
                
                hash_gt = {}
                
                try:
                    hash_gt['GT'] = r.genotype(sample)['GT']
                except:
                    pass

                try:
                    hash_gt['PL'] = ','.join(map(lambda p: str(p).replace('None','.'), r.genotype(sample)['PL']))
                except:
                    pass
                
                try:
                    AD = r.genotype(sample)['AD']
                    if type(AD) == int:
                        AD = (int(hash_gt['DP'])-AD,AD)
                    hash_gt['AD'] = ','.join(map(str,AD))
                except:
                    pass
                
                
                for ft in l_format:
                    ft_gt = hash_gt.get(ft,'.')
                    if ft == "GT" and str(ft_gt) == 'None':
                         ft_gt = "./."
                    elif ft == "AD" and (ft_gt == '.' or str(ft_gt) == 'None'):
                        ft_gt = "0,0"
                    elif ft == "PL" and (ft_gt == '.' or str(ft_gt) == 'None'):
                        ft_gt = "."
                    l_.append(ft_gt)
            
                gt_all = l_
                
                hash_variant = {}
                
                hash_variant['SNP_id'] = id
                hash_variant['qual'] = str(qual)
                hash_variant['filter'] = filter
                hash_variant['info'] = info
                hash_variant['format'] = ':'.join(l_format)
                hash_variant['format_info'] = ':'.join(gt_all)
                
                hash_table[(chr,pos,ref,alt)] = hash_variant
                
        except:
            
            raise RuntimeError('mosaic_somatic_variants_calling.parse_tissue_vcf: Some error has occurred in variant line %i:\n%d' % (i))    

    return hash_table,header,sample

########################################################################################################################
# functions that replaces in a given file, the pattern for the substitution
def replace(file_path, pattern, subst):
    #Create temp file
    fh, abs_path = mkstemp()
    with open(abs_path,'w') as new_file:
        with open(file_path) as old_file:
            for line in old_file:
                new_file.write(line.replace(pattern, subst))
    close(fh)
    #Remove original file
    remove(file_path)
    #Move new file
    move(abs_path, file_path)
    
########################################################################################################################
def __get_root_name(l_files):
    
    l_basename = map(lambda p: os.path.basename(p) , l_files)
    
    if len(l_files) == 1:
        return os.path.splitext(l_basename[0]) 
    
    l_ = []    
     
    for x in izip_longest(*l_basename,fillvalue=None):
        l_elem = list(set(list(x)))
        
        if len(l_elem) == 1:
            l_.append(l_elem[0])
        else:
            break
        
    root = ''.join(l_)
    
    if len(root) > 1:
        if root[-1] == '.':
            root = root[:-1]
            
    return root

########################################################################################################################
# bcftools and samtools, some INDELs cannot annotate properly
# example: ref: ccctcctcctcctcctc
# alt1: ccctcctcctcctcctcctc (insertion ctc) must be ref (-), alt (CTC)
# alt2: ccctcctcctcctcc (deletion tc) must be ref (TC), alt (-)
# example2: ref: TTT
# alt:T must be ref (TT) and alt (-)
# example3: ref: T
# alt: TTT must be ref (-), alt (TT) 


def check_minimal_INDEL_representation(pos,ref,alt):
    # insertion
    if ref.find(alt):
        diff = len(ref) - len(alt)
        diff = abs(diff) + 1
        alt_new = ref[-diff:]
        ref_new = alt_new[0]

    # deletion    
    elif alt.find(ref):
        diff = len(alt) - len(ref)
        diff = abs(diff) + 1
        ref_new = ref[-diff:]
        alt_new = ref_new[0]

    else:
        alt_new = alt
        ref_new = ref
        
    return pos,ref_new,alt_new


########################################################################################################################
# McArthur lab

def get_minimal_representation(pos, ref, alt): 
    # If it's a simple SNV, don't remap anything
    if len(ref) == 1 and len(alt) == 1: 
        return pos, ref, alt
    else:
        # strip off identical suffixes
        while(alt[-1] == ref[-1] and min(len(alt),len(ref)) > 1):
            alt = alt[:-1]
            ref = ref[:-1]
        # strip off identical prefixes and increment position
        while(alt[0] == ref[0] and min(len(alt),len(ref)) > 1):
            alt = alt[1:]
            ref = ref[1:]
            pos += 1
        #print 'returning: ', pos, ref, alt
        return pos, ref, alt 
    
    
########################################################################################################################
    
def __get_header_v2(vcf_file):
    
    l_header = []
    header   = ""
    
    fi = open(vcf_file,'r')
    
    for l in fi.readlines():
        
        if l[0] == '#':
            l_header.append(l)
        else:
            break
    
    fi.close()
    
    header = "".join(l_header)
    
    return header    

########################################################################################################################
def __get_header(record):
    
    header = ""
    
    try:
        hash_metadata = record.metadata
    except AttributeError:
        hash_metadata = {}
    
    for key in hash_metadata:
        if (type(hash_metadata[key])==str):
            header += "##%s=%s\n" % (key,hash_metadata[key])
        else:
            for met in hash_metadata[key]:
                if type(met) == str:
                    header += "##%s=%s\n" % (key,met)
                else:
                    header += "##%s=<" % (key)
                    for key_k in met: 
                        header += "%s=%s," % (key_k,str(met[key_k]))
                    header = header[:-1]
                    header += ">\n"
    
    #_Info = collections.namedtuple('Info', ['id', 'num', 'type', 'desc'])
    try:
        hash_info = record.infos
    except AttributeError:
        hash_info = {}
        
    for info in hash_info.keys():
        
        num_str = ""
        
        if hash_info[info].num == None:
            num_str = '.'
        else:
            num_str = str(hash_info[info].num)
        
        header += "##INFO=<ID=%s,Number=%s,Type=%s,Description=\"%s\">\n" % (hash_info[info].id,num_str,hash_info[info].type,hash_info[info].desc)
        
    #_Filter = collections.namedtuple('Filter', ['id', 'desc'])
    try:
        hash_filter = record.filters
    except AttributeError:
        hash_filter = {}
    
    for filt in hash_filter.keys():
        header += "##FILTER=<ID=%s,Description=\"%s\">\n" % (hash_filter[filt].id,hash_filter[filt].desc)
         
    #_Format = collections.namedtuple('Format', ['id', 'num', 'type', 'desc'])
    # It is set by hand 
    header += "##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">\n##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth (reads with MQ=255 or with bad mates are filtered)\">\n##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification\">\n"
    
    #'Contig', ['id', 'length']
    try:
        hash_contigs = record.contigs
    except AttributeError:
        hash_contigs = {}
    
    for cont in hash_contigs:
        header += "##contig=<ID=%s,length=%d>\n" % (hash_contigs[cont].id,hash_contigs[cont].length)  
    
    return header

########################################################################################################################
# Function that splits each multi-allelic row into 2 or more separates (DP and PL as well)
def split_multiallelic_vcf_to_simple(l_vcf,bcftools_path,threshold,logger):

    l_vcf_split = []
    
    logger.info('Splitting each variant with multi-allelic locus ...\n')
    
    for index, vcf_file in enumerate(l_vcf):

        # Extract the sample name from the VCF file
        bcftools_query_args = [bcftools_path, 'query', '-l', vcf_file]

        try:

            sample_name = subprocess.check_output(bcftools_query_args).rstrip()

        except subprocess.CalledProcessError, e:

            msg = 'variant_allele_fraction_filtering: Error running bcftools query -l :%s' % str(e)
            print msg
            raise RuntimeError(msg)

        fileName, fileExtension = os.path.splitext(vcf_file)
        vcf_file_split = fileName + '.split.vcf' 
        
        
        vcf_reader = vcf.Reader(filename=vcf_file)
         
        ## INFO fields
        hash_info = dict(vcf_reader.infos)
     
        ## header 
        header = __get_header(vcf_reader)        
         
        l_format = ['GT', 'PL', 'AD']
         
        fileName, fileExtension = os.path.splitext(vcf_file)
        vcf_file_split = fileName + '.split.vcf' 
         
        ## 1. Write in the corresponding splitted vcf file, the header
        fo = open(vcf_file_split,'w')
        fo.write(header)
 
        fo.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n" %sample_name)
         
        ## 2. For each vcf file (sample), iterate for each record, and separate in 2 o more rows, the multiallelic sites
         
        for i,r in enumerate(vcf_reader):
            if type(r) == _Record:
                 
                chr = r.CHROM
                pos = r.POS
                id = r.ID
                if id == None:
                    id = '.'
                ref = str(r.REF).strip()               
                alt = r.ALT
                 
                qual = r.QUAL                
                 
                info = r.INFO
                format = r.FORMAT               
                filter = r.FILTER
                 
                if filter == []:
                    filter = "PASS"
                else:
                    filter = str(filter[0])
                 
                l_s = r.samples # samples of each vcf
 
                # en este caso solo hay un sample
                for sample in l_s:
                    if 'GT' in format:
                        sample_GT = sample['GT']
                        if sample_GT is None:
                            continue
                    else:
                        sample_GT = '.'
                    if 'PL' in format:
                        l_PL = sample['PL']
                        #l_PL = sample_PL.split(",")
                    else:
                        sample_PL = '.'
                    if 'AD' in format:
                        sample_AD = sample['AD']
                    else:
                        sample_AD = '.'
                         
                l_alt = []
                if len(alt) > 1:                        
                    for x in alt:
                        l_alt.append(str(x))
                    alt_str = ','.join(l_alt)
                else:
                    alt_str = str(alt[0])
                     
                if len(l_alt) > 1:
                    
                    mavaf = info['MAVAF']
                    l_avaf = info['AVAF'].split('-')
                    ref_AD = sample_AD[0]
                     
                    for index,a in enumerate(l_alt):
                    
                        ## Split from the INFO field:AVAF. And if update the MAVAF value in the alternative allele cases
                        avaf = l_avaf[index]
                        mavaf = avaf
                         
                        alt_aux = str(alt[index])
                         
                        # l_PL extract the first 3 elements...
                        start = index * 3 
                        end = start + 3
                        l_PL_aux = l_PL[start:end]
                         
                        sample_AD_aux = [ref_AD,sample_AD[index+1]]
         
                        str_to_append = ""
                        try:
                            l_keys = list(set(info.keys()).difference(set(['AVAF','MAVAF'])))
                            for k in l_keys:
                                value = info[k]
                                if type(value) == list:
                                    value = ','.join(map(str,value))
                                str_to_append += ";" + "%s=%s" % (k,value)
                        except:
                            pass
                         
                         
                        new_info    = "AVAF=%s;MAVAF=%s" % (avaf,mavaf)
                         
                        info_print = new_info + str_to_append
                        if info_print[-1] == ";":
                            info_print = info_print[:-1]
                                 
                        # ACHTUNG!! si avaf < mosaicism, hay que cambiar la etiqueta de FILTER
                        if float(avaf) < threshold:
                            filter = "LowDepthAltMosaicism"
                         
                         
                        ## Build the FORMAT field
                        sample_PL = l_PL_aux
                        format_info = sample_GT + ":" + ','.join(map(str,sample_PL)) + ":" + ','.join(map(str,sample_AD_aux))
                         
                        # minimum INDEL representation
                        if r.INFO.has_key('INDEL'):                            
                            pos_new,ref_new,alt_new = check_minimal_INDEL_representation(pos,ref,a)
                            pos,ref,alt_aux = pos_new,ref_new,alt_new
 
                        alt_aux = alt_aux.strip()
                         
                        #CHROM    POS    ID    REF    ALT    QUAL    FILTER    INFO    FORMAT    SAMPLE_NAME
                        fo.write("%s\t%d\t%s\t%s\t%s\t%1.2f\t%s\t%s\t%s\t%s\n" % (chr,pos,id,ref,alt_aux,qual,filter,info_print,':'.join(l_format),format_info))
                     
                     
                else:
                     
                    # we do printing without doing anything: INFO and FORMAT fields tal y como estan
                    sample_PL = l_PL
                    format_info = sample_GT + ":" + ','.join(map(str,sample_PL)) + ":" + ','.join(map(str,sample_AD))
 
                    str_to_append = ""
                    try:
                        l_keys = list(set(info.keys()))
                        for k in l_keys:
                            value = info[k]
                            if type(value) == list:
                                value = ','.join(map(str,value))
                            str_to_append += ";" + "%s=%s" % (k,value)
                    except:
                        pass
                     
                    info_print = str_to_append
                    if info_print[-1] == ";":
                        info_print = info_print[:-1]
                    if info_print[0] == ";":
                        info_print = info_print[1:]
                     
                    # minimum INDEL representation
                    if r.INFO.has_key('INDEL'):
                        pos_new,ref_new,alt_new = check_minimal_INDEL_representation(pos,ref,alt_str)
                        pos,ref,alt_str = pos_new,ref_new,alt_new                    
                     
                    alt_str = alt_str.strip()
                    fo.write("%s\t%d\t%s\t%s\t%s\t%1.2f\t%s\t%s\t%s\t%s\n" % (chr, pos, id, ref, alt_str, float(qual), filter, info_print, ':'.join(l_format),format_info))
                                 
        fo.close()
         
        logger.info("%s split done \n" %(vcf_file_split))    
        
        l_vcf_split.append(vcf_file_split)

    logger.info("Splitting finished! ...\n")
    
    return l_vcf_split

########################################################################################################################
# This function computes the AVAF and checks whether a record's VAF is >= threshold (0.07 by default)
# min_allele_balance: the minimum number account in the alternative allele

def variant_allele_fraction_filtering(l_vcf,bcftools_path,threshold,logger,min_allele_balance=0.02):

    l_vcf_avaf_filtered = []
    logger.info('AVAF new INFO field determination...\n')
    
    for index, vcf_file in enumerate(l_vcf):

        # Extract the sample name from the VCF file
        bcftools_query_args = [bcftools_path, 'query', '-l', vcf_file]

        try:

            sample_name = subprocess.check_output(bcftools_query_args).rstrip()

        except subprocess.CalledProcessError, e:

            msg = 'variant_allele_fraction_filtering: Error running bcftools query -l :%s' % str(e)
            print msg
            raise RuntimeError(msg)

        fileName, fileExtension = os.path.splitext(vcf_file)
        vcf_file_filtered = fileName + ".AVAF.filtered" + str(int(threshold*100)) + "%.vcf"  # vcf file with the AVAF filtered depending on the threshold 
        
        dict_filters = {}
        
        vcf_reader = vcf.Reader(filename=vcf_file)
        dict_filters.update(vcf_reader.filters)
        
        ## 1. Define the new INFO field that will be included
        
        vcf_reader.infos['AVAF'] = _Info('AVAF', 1, 'String', 'Alternative variant allele frequency', '', '')
        vcf_reader.infos['MAVAF']  = _Info('MAVAF', 1, 'String', 'Maximum alternative variant allele frequency, the most representative one', '', '')
        #vcf_reader.infos['ALTSTR']  = _Info('ALTSTR',1,'String',"All the alternative alleles with a frequency higher than mosaicism threshold")

        vcf_reader.filters['LowDepthAltMosaicism'] = _Filter('LowDepthAltMosaicism', 'The most representative alternative variant allele frequency is lower than the mosaicism threshold')
        
        ## header 
        header = __get_header(vcf_reader)        
        
        l_format  = ['GT','PL','AD'] 
                
        
        ## 2. Write the header into the corresponding filtered VCF file
        fo = open(vcf_file_filtered, 'w')
        fo.write(header)
        
        fo.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n" % sample_name)
        
        # to break (continue) in inner for loop
        keeplooping = False
        
        ## 3. For each vcf file (sample), iterate for each record, and compute the AVAF (alternative variant allele fraction)
        for i,r in enumerate(vcf_reader):
            
            if type(r) ==  _Record:
                
                chr = r.CHROM
                pos = r.POS
                id = r.ID
                if id == None:
                    id = '.'
                ref = r.REF
                alt = r.ALT
                qual = r.QUAL
                
#                pos_new,ref_new,alt_new = get_minimal_representation(pos, ref, alt)
                
#                pos,ref,alt = pos_new,ref_new,alt_new     
                
                # si hay presencia de un alelo alternativo seguimos, si no, lo quitamos del vcf (no hay cambio), limpio
                # alt = [None]
                
                if alt != [None]:
                    l_alt = []
                    if len(alt) > 1:                        
                        for x in alt:
                            l_alt.append(str(x))
                        alt_str = ','.join(l_alt)
                    else:
                        alt_str = str(alt[0])
                        l_alt = alt_str
                                        
                    info = r.INFO
                    format = r.FORMAT
                    
                    l_s = r.samples # samples of each vcf
                    
                    # en este caso solo hay un sample
                    for sample in l_s:
                        if 'GT' in format:
                            sample_GT = sample['GT']
                            if sample_GT is None:
                                keeplooping = True
                            if sample_GT == "./.":
                                keeplooping = True
                        else:
                            sample_GT = '.'
                        if 'PL' in format:
                            l_PL = sample['PL']
                            #l_PL = sample_PL.split(",")
                        else:
                            sample_PL = '.'
                        if 'AD' in format:
                            sample_AD = sample['AD']
                        else:
                            sample_AD = '.'
                    
                    if keeplooping:
                        continue
                    
                    ref_depth = float(sample_AD[0])
                    if len(sample_AD) > 2:
                        l_alt_depth = sample_AD[1:len(sample_AD)]
                        alt_depth = l_alt_depth[0]
                    else:
                        l_alt_depth = float(sample_AD[1])
                        alt_depth = sample_AD[1]
                        
                    
                    # If the site is homo alt ==> ref_depth could be 0
                    alt_major = False
                    
                    if ref_depth == 0:
                        
                        ref_depth = alt_depth
                             
                    elif alt_depth > ref_depth:
                        
                        alt_major = True
                    
                    if len(sample_AD) > 2:
                        
                        if len(sample_AD) == 3:
                            alt_depth2 = float(l_alt_depth[1])
                            
                            total_depth = alt_depth + ref_depth
                            avaf = float(alt_depth/total_depth)
                            
                            if alt_major :
                                total_depth2 = alt_depth2 + alt_depth
                            else:
                                total_depth2 = alt_depth2 + ref_depth
                                
                            avaf2 = float(alt_depth2/total_depth2)
                            
                            # filtramos alt1 ? si alt1 se filtra tambien el alt2 (proque van en orden)
                            if avaf <= min_allele_balance:
                                continue
                            
                            # si alt2 se filtra, quitamos del alt_str
                            if avaf2 <= min_allele_balance:
                                l_alt.pop()
                                sample_AD.pop()
                                # eliminamos los 3 ultimos GT's correspondientes al alelo     
                                del l_PL[-3:]                        
                                alt_str = ','.join(l_alt) 
                                
                        
                        elif len(sample_AD) == 4:
                            alt_depth2 = float(l_alt_depth[1])
                            alt_depth3 = float(l_alt_depth[2])
                            
                            total_depth = alt_depth + ref_depth
                            avaf = float(alt_depth/total_depth)
                            
                            if alt_major:
                                total_depth2 = alt_depth2 + alt_depth
                            else:
                                total_depth2 = alt_depth2 + ref_depth
                                
                            avaf2 = float(alt_depth2/total_depth2)
                            
                            if alt_major:
                                total_depth3 = alt_depth3 + alt_depth
                            else:
                                total_depth3 = alt_depth3 + ref_depth
                                
                            avaf3 = float(alt_depth3/total_depth3)
                            
                            # filtramos alt1??? si alt1 se filtra tambien el alt2 y alt3, todo.
                            if avaf <= min_allele_balance:
                                continue

                            # si alt2 se filtra, quitamos del alt_str
                            if avaf2 <= min_allele_balance:
                                # solo seria al primer alt, porque avaf3 sera menor (van en orden)
                                l_alt = l_alt[0]
                                sample_AD = sample_AD[0:2]
                                l_PL = l_PL[0:3]
                                alt_str = str(l_alt) 
                                
                            elif avaf3 <= min_allele_balance:
                                l_alt.pop()
                                sample_AD.pop()
                                del l_PL[-3:]
                                alt_str = ','.join(l_alt)
                            
                            
                    else:    
                        total_depth = alt_depth + ref_depth
                        avaf = float(alt_depth/total_depth)
                        
                        # si avaf < min_allele_balance directamente lo filtramos. No pasa ni a la etiqueta de lowMosaicism. kaka
                        if avaf <= min_allele_balance:
                            continue
    
                    ## l_PL must be a 3 multiple ==> 3 * len(l_alt). We must remove or add 255 (for instance, we do not use this data) not to have errors during the annotation (otherwise, snpEff crashes)
                    if not isinstance(l_alt,str):
                        if len(l_PL) != 3*len(l_alt):
                            diff = len(l_PL) - (3 * len(l_alt))
                            if diff > 0:
                                del l_PL[-diff:]
                            else:
                                l_PL.extend(np.repeat(255, np.abs(diff)))
                    else:
                        if len(l_PL) != 3:
                            diff = len(l_PL) - 3
                            if diff > 0:
                                del l_PL[-diff:]
                            else:
                                l_PL.extend(np.repeat(255, np.abs(diff)))
                    
                    ## MAVAF
                    if len(sample_AD) > 2:
                        # we identify the most representative alternative allele, and we look whether is >= threshold mosaicism
                        max_allele = float(max(map(int,l_alt_depth))) #  a.index(max(a))
                        mavaf = float(max_allele/(ref_depth + max_allele))
                        
                        if mavaf >= threshold:
                            filter = "PASS"
                        else:
                            filter = "LowDepthAltMosaicism"
                        
                    else:
    
                        mavaf = avaf
                        if mavaf >= threshold:
                            filter = "PASS"
                        else:
                            #   ##FILTER=<ID=LowDepthAltMosaicism,Description="Allele depth lower than the mosaic threshold">
                            filter = "LowDepthAltMosaicism"

                    
                    ##AVAF and malt. We now define malt, maximum alternative allele (the one with the major frequency).
                    ##alt tendra los alelos alternativos cuya frecuencia es > threshold mosaicism
                    
                    if type(l_alt) == list :
                        if len(l_alt) == 3:
                            avaf_str = str(round(avaf,2)) + "-" + str(round(avaf2,2)) + '-' + str(round(avaf3,2))
                            malt = alt_str[0]
    
                            if len(l_alt) > 0:
                                malt = l_alt[0]
                                alt = ','.join(l_alt)
                                
                            else:
                                malt = str(alt[0])
                                alt = str(l_alt[0])
                                
                        elif len(l_alt) == 2:
                            avaf_str = str(round(avaf,2)) + "-" + str(round(avaf2,2))
                            
                            if len(l_alt) > 0:
                                malt = l_alt[0]
                                alt = ','.join(l_alt)
                            
                            else:
                                malt = str(alt[0])
                                alt = str(alt[0])

                        else:
                            avaf_str = str(round(avaf,2))
                            
                            if type(l_alt) == list:
                                malt = l_alt[0]
                                alt = ','.join(l_alt)
                            else:
                                malt = str(alt[0])
                                alt = str(alt[0])

                    
                    else:
                        avaf_str = str(round(avaf,2))
                        
                        if type(l_alt) == list:
                            malt = l_alt[0]
                            alt = ','.join(l_alt)
                        else:
                            malt = str(alt[0])
                            alt = str(alt[0])


                    
                    ## INFO new fields generation
                    new_info    = "AVAF=%s;MAVAF=%s" % (avaf_str,str(round(mavaf,2)))
                    #new_info    = "AVAF=%s;MAVAF=%s;ALTSTR=%s" % (avaf_str,str(round(mavaf,2)),alt_str)
    
                    str_to_append = ""
                    try:
                        l_keys = list(set(info.keys()).difference(set(['AVAF','MAVAF'])))
                        for k in l_keys:
                            value = info[k]
                            if type(value) == list:
                                value = ','.join(map(str,value))
                            str_to_append += ";" + "%s=%s" % (k,value)
                    except:
                        pass
                    
                    info_print = new_info + str_to_append
                    if info_print[-1] == ";":
                        info_print = info_print[:-1]
                        
                    
                        
                    ## Build the FORMAT field
                    #sample_PL = ','.join(map(str,l_PL))
                    sample_PL = l_PL
                    format_info = sample_GT + ":" + ','.join(map(str,sample_PL)) + ":" + ','.join(map(str,sample_AD))
                    
            
                    #CHROM    POS    ID    REF    ALT    QUAL    FILTER    INFO    FORMAT    SAMPLE_NAME
                    fo.write("%s\t%d\t%s\t%s\t%s\t%1.2f\t%s\t%s\t%s\t%s\n" % (chr,pos,id,ref,alt,qual,filter,info_print,':'.join(l_format),format_info))
        fo.close()
        logger.info("%s created \n" %(vcf_file_filtered))
        
        l_vcf_avaf_filtered.append(vcf_file_filtered)
    
    return l_vcf_avaf_filtered

########################################################################################################################

def capture_callings(l_bcf,bcftools_path,logger):
    
    l_vcf = []
    
    for bcf_file in l_bcf:
        
        fileName, fileExtension = os.path.splitext(bcf_file)
        vcf_file = fileName + ".vcf"
        
        logger.info("bcftools call --> %s\n" %(bcf_file))
        #$bcf_path call -Am -O v $bcf_file -o $vcf_file
        bcftools_args = [bcftools_path,'call','-Am','-O','v',bcf_file,'-o',vcf_file]
        bcftools_sal = Popen(bcftools_args,stdin=PIPE,stdout=PIPE,stderr=PIPE,close_fds=True,bufsize=1)
        (trash,logdata) = bcftools_sal.communicate()
        bcftools_sal.wait()
        
        if logdata != "":
            if logdata.lower().find("error") != -1:
                raise RuntimeError('mosaic_somatic_variants_calling.capture_callings: Error in bcftools call:\n%s\n' % (logdata))
    
        logger.info("bcftools calling done!\n")
        
        l_vcf.append(vcf_file)
        
        
    return l_vcf


########################################################################################################################
def check_samtools_bctools_versions(samtools_path,bcftools_path):
    
    # /home/kibanez/Escritorio/samtools-1.3.1/samtools --version    
    samtools_args = [samtools_path,'--version']
    file_s = "samtools_version"
    f_samtools = open(file_s,'w')
    samtools_sal = Popen(samtools_args,stdin=PIPE,stdout=f_samtools,stderr=PIPE,close_fds=True,bufsize=1)
    (trash,logdata) = samtools_sal.communicate()
    samtools_sal.wait()
    
    if logdata != "":
        if logdata.lower().find("error") != -1:
            raise RuntimeError('mosaic_somatic_variants_calling.check_samtools_bcftools_versions: Error in samtools --version:\n%s\n' % (logdata))
    
    f_samtools.close()
    
    # /home/kibanez/Escritorio/bcftools-1.3.1/bcftools --version
    bcftools_args = [bcftools_path,'--version']
    file_b = "bcftools_version"
    f_bcftools = open(file_b,'w')
    bcftools_sal = Popen(bcftools_args,stdin=PIPE,stdout=f_bcftools,stderr=PIPE,close_fds=True,bufsize=1)
    (trash,logdata) = bcftools_sal.communicate()
    bcftools_sal.wait()
    
    if logdata != "":
        if logdata.lower().find("error") != -1:
            raise RuntimeError('mosaic_somatic_variants_calling.check_samtools_bcftools_versions: Error in bcftools --version:\n%s\n' % (logdata))
    
    f_bcftools.close()
    
    # We compare the versions' output
    v_samtools = ""
    with open(file_s) as infile:
        for line in infile:
            if 'samtools' in line:
                v_samtools = line
    
    v_bcftools = ""
    with open(file_b) as infile:
        for line in infile:
            if 'bcftools' in line:
                v_bcftools = line
    
    
    if v_samtools.split(" ")[1] == v_bcftools.split(" ")[1]:
        versions = True
    else:
        versions = False
    
    
    # we do remove the temporal files
    if os.path.isfile(file_s):
        os.remove(file_s)
    if os.path.isfile(file_b):
        os.remove(file_b)
        
    return versions



########################################################################################################################

def mpileup_calling(l_bam,bed,fasta,samtools_path,variant_path,q_value,q_value2,logger):

    l_bcf = []
    
    for i, bam_file in enumerate(l_bam):
        
        logger.info('Samtools mpileup --> %s \n' %bam_file)

        bam_name = os.path.splitext(os.path.basename(bam_file))[0]

        sample_folder = os.path.join(variant_path, bam_name)
        if not os.path.exists(sample_folder):
            os.mkdir(sample_folder)

        bcf_file = os.path.join(sample_folder, bam_name + '.bcf')

        f_bcf = open(bcf_file,'w')
        
        #$samtools_path mpileup -Buf $ref_fasta --positions $bed --output-tags DP,AD,ADF,ADR,SP,INFO/AD,INFO/ADF,INFO/ADR -q 0 -Q 0 $bam_file > $bcf_file
        mpileup_args = [samtools_path, 'mpileup', '-Buf', fasta, '--positions', bed, '--output-tags', 'DP,AD,ADF,ADR,SP,INFO/AD,INFO/ADF,INFO/ADR',
                        '-q', q_value, '-Q', q_value2, bam_file]
        mpileup_sal = Popen(mpileup_args, stdin=PIPE, stdout=f_bcf, stderr=PIPE, close_fds=True, bufsize=1)
        (trash, logdata) = mpileup_sal.communicate()
        mpileup_sal.wait()
        
        if logdata != "":
            if logdata.lower().find("error") != -1:
                raise RuntimeError('mosaic_somatic_variants_calling.mpileup_calling: Error in samtools mpileup:\n%s\n' % logdata)
        
        f_bcf.close()
        
        l_bcf.append(bcf_file)
    
        logger.info('samtools mpileup done!\n')
        
    return l_bcf


########################################################################################################################
# Function that selectes the INFO and FORMAT fields, and then, downgrades to VCF4.1. And it also changes the FORMAT AD field, that is badly defined in the new version VCF4.2
def select_columns_and_vTransform(l_vcf,bcftools_path,logger):
    l_vcf_filtered = []
    logger.info("Extracting fields from the FORMAT and INFO fields and downgrading to VCF4.1 ...\n")
    
    for vcf_file in l_vcf:
        fileName, fileExtension = os.path.splitext(vcf_file)
        vcf_filtered = fileName + "_filtered.vcf"
        
                    
        # 1) Extract FORMAT and INFO fields        
        # $bcf_path annotate -x INFO/ADF,INFO/ADR,INFO/SGB,INFO/AD,INFO/RPB,INFO/MQB,INFO/MQSB,INFO/BQB,INFO/MQ0F,INFO/AC,INFO/AN,INFO/DP4,INFO/VDB,^FORMAT/GT,FORMAT/AD,FORMAT/PL $vcf_file -o $vcf_filtered
        
        manipulate_sal = Popen([bcftools_path,'annotate','-x','INFO/ADF,INFO/ADR,INFO/SGB,INFO/AD,INFO/RPB,INFO/MQB,INFO/MQSB,INFO/BQB,INFO/MQ0F,INFO/AC,INFO/AN,INFO/DP4,INFO/VDB,^FORMAT/GT,FORMAT/AD,FORMAT/PL',vcf_file,'-o',vcf_filtered],stdin=PIPE, stdout=PIPE, stderr=PIPE,close_fds=True,bufsize=1)
        (trash,logdata) = manipulate_sal.communicate()
        manipulate_sal.wait()
        
        if logdata != "":
            if logdata.lower().find('fallo') != -1:
                raise RuntimeError('mosaic_somatic_variants_calling.clean_columns_and_vTransform: Error when launching bcftools annotate -x :\n%s' % (logdata))
        
        
        # 2) Downgrade VCF4.2 version to VCF4.1
        # sed "s/##fileformat=VCFv4.2/##fileformat=VCFv4.1/" $vcf_filtered | sed "s/(//" | sed "s/)//" | sed "s/,Version=\"3\">/>/" | sed 's/##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">/##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed"> /g' > $vcf_filtered_v4_1
        
        if not os.path.isfile(vcf_filtered):
            
            raise IOError('The vcf file already filtered does not exist %s' % (vcf_filtered))
        
        else:        
        
            replace(vcf_filtered,"##fileformat=VCFv4.2","##fileformat=VCFv4.1")
            replace(vcf_filtered,"(","")
            replace(vcf_filtered,")","")
            replace(vcf_filtered,",Version=\"3\">",">")
            replace(vcf_filtered,"##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">","##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">")
        
        l_vcf_filtered.append(vcf_filtered)
    
    logger.info("done!\n")
            
    return l_vcf_filtered

########################################################################################################################

def run(argv=None):
    
    if argv is None: argv = sys.argv    
   
    parser = OptionParser(add_help_option=True,description="",formatter= optparse.TitledHelpFormatter(width = 200))
    
    parser.add_option("--cfg", default=None,
                      help= 'Input cfg file',
                      dest="f_cfg")
    parser.add_option("--m", default=None,
                      help= 'The desired threshold of the mosaicism (by default the minimum threshold is set to 0.07 (7%))',
                      dest="f_threshold")
    parser.add_option("--s", default=False, action="store_true",
                      help= 'The cfg file must contain tissue, blood if it is available and control samples in order to analyze the mosaicism somatic changes',
                      dest="f_somatic")

                        
    (options, args) = parser.parse_args(argv[1:])

    if len(argv) == 1:
        sys.exit(0)
    
    parser.check_required("--cfg")
    
    formatter = logging.Formatter('%(asctime)s - %(module)s - %(levelname)s - %(message)s')
    console = logging.StreamHandler()
    console.setFormatter(formatter)
    console.setLevel(logging.INFO)
    logger = logging.getLogger("preprocess")
    logger.setLevel(logging.INFO)
    logger.addHandler(console)

    try:
        cfg_file = options.f_cfg
         
        if not os.path.exists(cfg_file):
            raise IOError('NGS_call_somatic_mosaic_variants: The cfg file %s does not exist' % cfg_file)
        
        hash_cfg = mod_cfg.read_cfg_file(cfg_file)
        
        f_threshold = options.f_threshold
        
        if f_threshold is None:
            
            f_threshold = 0.07
            
        else:
            
            f_threshold = float(f_threshold)
        
        # input required
        analysis_bed = hash_cfg.get('analysis_bed', '')
        input_files = hash_cfg.get('input_files', '')
        q_value = hash_cfg.get('q_value', '20')
        q_value2 = hash_cfg.get('Q_value', '20')
        
        # reference required
        ref_fasta = hash_cfg.get('ref_fasta', '')

        # software required
        samtools_path = hash_cfg.get('samtools_path', '')
        bcftools_path = hash_cfg.get('bcftools_path', '')
        
        # output required
        variant_path = hash_cfg.get('variant_path', '')
        
        
        if q_value is None:
            logger.info('The value of the base threshold quality (q) would be 0 by default. Please, introduce a proper q_value \n')

        if q_value2 is None:
            logger.info('The value of the mapping quality threshold (Q) would be 0 by default. Please, introduce a proper Q_value \n')
        
        if not os.path.isfile(ref_fasta):
            raise IOError("The file does not exist. %s" % ref_fasta)

        if not os.path.isfile(analysis_bed):
            raise IOError("The analysis bed does not exist. %s" % analysis_bed)
        
        if not os.path.exists(samtools_path):
            raise IOError('The samtools_path path does not exist %s' % samtools_path)

        if not os.path.exists(variant_path):
            raise IOError('The variant path does not exist %s' % variant_path)
        

        f_somatic = options.f_somatic
        
        if f_somatic == False:
            
            # if '--s' option is not given, then the mosaic variants are computed from the aligned bam files
            l_bam = input_files.split(",")
            
            for i in l_bam:
                if not os.path.isfile(i):
                    raise IOError("NGS_call_somatic_mosaic_variants: The given BAM file does not exist. %s" % i)
            
            
            # IMPORTANT PRELIMINAR STEP: check whether the samtools and bcftools versions are the same (it usually updates at the same level and time)
            versions_tools = check_samtools_bctools_versions(samtools_path,bcftools_path)
            
            if not versions_tools:
                raise IOError('Warning: samtools and bcftools versions must be the same. It is absolutely required for an adequate usage of this tool. Please visit http://www.htslib.org/download/ \n')
    
            logger.info("The mosaicism variant detection starts...\n")
            
            # (1) mpileup the BAM files
            l_bcf = mpileup_calling(l_bam, analysis_bed, ref_fasta, samtools_path, variant_path, q_value, q_value2, logger)
             
            # (2) capture of the calls        
            l_vcf = capture_callings(l_bcf, bcftools_path, logger)
             
            # (3) extract the desired columns and transform them from VCF4.1 to VCF4.2 in order to work easily with them
            l_vcf_trans = select_columns_and_vTransform(l_vcf, bcftools_path, logger)
             
            # (4) manual filtering: creation of AVAF (alternative variant allele fraction) depending on the threshold (--m) inserted
            l_vcf_avaf = variant_allele_fraction_filtering(l_vcf_trans, bcftools_path, f_threshold, logger)
             
            # (5) samples_run/controls annotation
             
            # (5.1) split multi-allelic sites into different simple rows
            l_vcf_split = split_multiallelic_vcf_to_simple(l_vcf_avaf, bcftools_path, f_threshold, logger)
                 
            # (5.2) combine: samples_run (compute how many times the variant is detected among the samples included in the cfg file
            mod_variant.annotate_vcf(l_vcf_split, variant_path, hash_cfg)
            
            logger.info('The mosaicism variant detection finished!\n')
            logger.info('You only need to annotate the resulting VCF files ;_)\n')
    
        else:
            
            logger.info('The somatic mosaicism variant detection starts...\n')
            
            # input files
            tissue_id = hash_cfg.get('tissue_id', '')
            blood_id = hash_cfg.get('blood_id', '')
            control_id = hash_cfg.get('control_id', '')
            sample_name = hash_cfg.get('sample_name', '')
            
            if not os.path.isfile(tissue_id):
                raise IOError("The variant file corresponding to the tissue does not exist. %s" % tissue_id)
            
            if blood_id != "":
                if not os.path.isfile(blood_id):
                    raise IOError("The variant file corresponding to the blood does not exist. %s" % blood_id)
                
            elif blood_id is None:
                warnings.warn('The variant file corresponding to the blood is not given. Thus, the somatic analysis would be done filtering with the given control samples. \n')
                
            l_control = control_id.split(",")
            for i in l_control:
                if not os.path.isfile(i):
                    raise IOError("The variant file corresponding to on of the control does not exist. %s" % i)

            if sample_name == "":
                raise IOError("The sample name must be have info. %s" % i)
            
            fileName, fileExtension = os.path.splitext(tissue_id)
            file_somatic = fileName + "_VS_blood_and_controls"
                     
            hash_tissue,header,sample = parse_tissue_vcf(tissue_id)
            
            if hash_tissue == {}:
                return hash_tissue
            
            # Sometimes there is no paired-blood sample which corresponds to the tissue sample. In those cases, we only filtered by the control samples
            if blood_id != "":
                
                hash_table_tmp = parse_control(blood_id, hash_tissue, blood=True)
            
                hash_tissue = hash_table_tmp

            for vcf_control in l_control:
                
                hash_table_tmp = parse_control(vcf_control, hash_tissue)
            
                hash_tissue = hash_table_tmp
            
            vcf_file_somatic = file_somatic + ".vcf"
            
            fo = open(vcf_file_somatic, 'w')
            fo.write(header)
            fo.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n" % sample)
            
            chr_cyto_sorted = map(lambda i: "chr%i" % i, range(1, 23))+['chrX', 'chrY', 'chrM']
            
            l_fields = ['qual', 'filter', 'info', 'format', 'format_info']
            
            for key in sorted(sorted(hash_tissue.keys(),key=itemgetter(1)),key=lambda x: chr_cyto_sorted.index(x[0])):        
                l_key = map(str,list(key))
                l_key.insert(2,hash_tissue[key].get('SNP_id','.'))
                key_str = '\t'.join(l_key)
                
                fo.write(key_str + '\t' + '\t'.join(map(lambda field: hash_tissue[key].get(field,'.'), l_fields))+'\n')
        
            fo.close()
            
            logger.info('The somatic mosaicism variant detection finished!\n')
            logger.info('You only need to annotate the resulting VCF file ;)\n')
            
            
    except:
        print >> sys.stderr, '\n%s\t%s' % (sys.exc_info()[0],sys.exc_info()[1])
        sys.exit(2)

########################################################################################################################
     
if __name__=='__main__':
    
    run()

