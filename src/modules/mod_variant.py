#!/usr/bin/python

import sys, re, os, string, urllib, time, math, random, subprocess, shutil, itertools

from os import path as osp

from os import stat

import optparse

from subprocess import Popen , PIPE

from itertools import izip_longest,chain,product

from operator import itemgetter

import itertools

import vcf

from vcf.parser import _Filter,_Record,_Info,_Format,_SV

from vcf.utils import walk_together

import pysam

import numpy,scipy

from collections import namedtuple 
 

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

# Function to split and include within a list, vcf's pass and filtered in each variant
def clean_walker(l_set_walker):
    l_clean = []
    
    for i in l_set_walker:
        if "set" in i:
            aux = i.split('=')
            l_clean.append(aux[1])
        else:
            l_clean.append(i)
    
    return l_clean
        
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
        
    if len(l_rod_priority) <> len(l_vcf_files):
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
    
    if mode <> "combine":

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
            
            if type(r_i) <> vcf.model._Record:
                continue
            
            chrm  = r_i.CHROM
            pos   = r_i.POS
            ref   = r_i.REF
            l_alt_r = map(lambda sb: sb.__repr__() , r_i.ALT)
            #l_s     = map(lambda s: s.sample, r_i.samples)
            
            l_s = []
            for c in r_i.samples:
                if c.gt_bases <> None:
                    l_s.append(c.sample)
            
            if l_alt_r == ['None']: # When the site has being for all the samples genotyped and is HOM_REF
                l_alt_r   = ['.']
                 
            l_alt.extend(l_alt_r)
            
            ### Build the ALT fields
            l_alt = list(set(l_alt))
            alt_str = ','.join(l_alt)
            if l_alt <> ['.']:
                if set(l_alt).intersection(set(['.']))<> set([]):
                #if filter(lambda a: a=='.', l_alt) <> []:
                    alt_str = alt_str.replace('.', '')
                    alt_str = alt_str.replace(',', '') 
            
                            
            hash_info = dict(r_i.INFO)
            
            is_HOM_REF = False
            is_MULTIALLELIC = False
            
            for sample in l_s:
                if r_i.genotype(sample)['GT'] == "0/0":
                    is_HOM_REF = True
                    num_calls += 1
                elif r_i.genotype(sample)['GT'] <> "0/0" and r_i.genotype(sample)['GT'] <> None:
                    num_calls    += 1
                    num_variants += 1
                     
            if len(l_alt_r) > 1:
                is_MULTIALLELIC = True
            
            if r_i.ID <> None:
                l_id.extend(r_i.ID.split(';'))
            else:
                if not is_HOM_REF:
                    l_id.append('.')
                elif is_MULTIALLELIC:
                    l_id.append('.')
                else:
                    l_id = ['.']
            
            if r_i.QUAL <> None:
                sum_QUAL += r_i.QUAL
            
            if r_i.QUAL <> None:
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
                
                if r_i.FILTER <> [] and r_i.FILTER <> None:
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
                if s<>None:
                    set_str += s+'--'
            set_str = set_str[:-2]
        
                
        ### Generating the string of the FILTER field. If any of the records with call is filtered, it is included in the FILTER string
        ### PASS means all the calls do PASS
        num_filtered = len(hash_filter.keys())
        l_filter = list(set(l_filter))
        filter_string  = "PASS"
        if l_filter <> []: # in case any call is filtered
            #if num_calls == num_filtered:
            # JUNE'16: OJO!!! if all the samples in the pool are filtered, so, filteredInALL and NC = num_samples. Otherwise, we do not know which samples are filtered
            if num_samples == num_filtered:
                set_str = "FilteredInAll"
                
            filter_string = ""
            for f in l_filter:
                if f<>'PASS':
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
                                
                            if filter_s <> 'PASS':
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
                           
                            if filter_s <> 'PASS' and l_set_walker == ["FilteredInAll"]:
                                l_set_walker = ["FilteredInAll"]
                            elif filter_s <> 'PASS' and l_set_walker <> ["FilteredInAll"]:
                                l_set_walker.append(l_set_unique)    
                            else:
                                l_set_walker = [l_set_unique] + l_set_walker
            
                set_str = '--'.join(l_set_walker)
                
                if set_str <> "FilteredInAll":
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
            if filter(lambda gt: gt == '0/0' or gt == 'None', l_gt) <> []:
                continue
        
        gt_all  = map(lambda i: ':'.join(i) , l_gt)
        ## Build the ID field
        l_id = filter(lambda id: id<>'.', list(set(l_id)))
        if l_id == []:
            l_id = ['.']
        id = ','.join(list(set(l_id)))
        
#         ### Build the ALT fields
#         l_alt = list(set(l_alt))
#         alt_str = ','.join(l_alt)
#         if l_alt <> ['.']:
#             if set(l_alt).intersection(set(['.']))<> set([]):
#             #if filter(lambda a: a=='.', l_alt) <> []:
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


    
def __change_sample_name(vcf_file,new_label='overlapping'):
    
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


def annotate_vcf(l_vcf_files_run,output_path,hash_cfg):
    
    """
        The goal of this function is to include the info of:
         - the variants of the same run (in order to detect possible artefacts) 
         - the variants of a set of controls (to disambiguate polymophisms from rare variants)
    """
    
    if hash_cfg.has_key('ref_fasta'):
        ref_fasta_file = hash_cfg['ref_fasta']
    else:
        raise AttributeError("mod_variant.variant_calling: The fasta file of the reference has not been configured properly") 
        
    if hash_cfg.has_key('gatk_path'):
        gatk_path = hash_cfg['gatk_path']
    else:
        raise AttributeError("mod_variant.variant_calling: The path of the gatk jar file has not been configured properly")
    
    l_vcf_files_run_out = []
    
    t1 = time.time()
    sys.stdout.write('Creating enriched annotations in the vcf files ...\n')
    sys.stdout.flush()
        
    #### 1. Get the variants that overlap in all the samples in a same run. The overlapping changes are written in the file 'vcf_overlapping_variants_within_run'
    sys.stdout.write('\nGetting the variants that completely overlap in all the samples of the run...\n')
      
    ###vcf_overlapping_variants_within_run = __annotate_overlapping_variants(l_vcf_files_run,output_path,hash_cfg)
    vcf_overlapping_variants_within_run = combine_vcf_v2(l_vcf_files_run,[],"combine",output_path)    
    
    ### the samples of the file 'vcf_overlapping_variants_within_run' are replaced to generic sample 'overlapping'
    vcf_tmp = __change_sample_name(vcf_overlapping_variants_within_run)
    
    name_vcf_overlapping_variants_within_run = os.path.join(output_path,'combined_overlapping_variants_within_run.vcf')
    if os.path.exists(vcf_tmp):
        os.rename(vcf_tmp,name_vcf_overlapping_variants_within_run)
         
    if os.path.exists(vcf_tmp+'.idx'):
        os.remove(vcf_tmp+'.idx')
         
    sys.stdout.write('[done]\n')
    sys.stdout.flush()
    
    ## 2. Combine the variants of controls
    sys.stdout.write('\nIncluding the information of controls...\n')
    sys.stdout.flush()
    
    l_vcf_files_controls = map(itemgetter(1),hash_cfg.get('controls',[]))
    
    if l_vcf_files_controls <> []:
        
        num_controls = len(l_vcf_files_controls) 
        
        sys.stdout.write('Number of controls loaded: %s\n' % (num_controls))
        
        l_controls = map(lambda i: 'control_%d' % (i), range(num_controls))
    
        controls_vcf_old = combine_vcf_v2(l_vcf_files_controls,l_controls,"combine",output_path)
        ###controls_vcf_old = combine_vcf(l_vcf_files_controls,l_controls,ref_fasta_file,gatk_path,"multiple",output_path)
    
        controls_vcf = os.path.join(output_path,'combined_controls.vcf')
    
        
        if os.path.exists(controls_vcf_old):
            os.rename(controls_vcf_old,controls_vcf)
         
        if os.path.exists(controls_vcf_old+'.idx'):
            os.remove(controls_vcf_old+'.idx')
            
        control_reader = vcf.Reader(filename=controls_vcf)
        run_reader = vcf.Reader(filename=name_vcf_overlapping_variants_within_run)
        
        # march'16: we do now hash with 4-tuple id's: chrom, pos, ref and alt (because combine 'combine' puts in 2 different lines identical chr-pos but with different ref-alt info!!
        #hash_control = dict(map(lambda r: ((r.CHROM,r.POS,r.REF),r), control_reader))
        
        hash_control     =      map( lambda r: (
                                            map ( lambda a: ((r.CHROM,r.POS,r.REF,str(a)),r),
                                                  r.ALT), 
                             )[0],
                             control_reader
                            
                           )
        hash_control2 = list(itertools.chain.from_iterable(hash_control))
        
        hash_control = dict(hash_control2)
        
                
        hash_run     =      map( lambda r: (
                                            map ( lambda a: ((r.CHROM,r.POS,r.REF,str(a)),r),
                                                  r.ALT), 
                             )[0],
                             run_reader
                            
                           )
                        
        hash_run2 = list(itertools.chain.from_iterable(hash_run))
        
        hash_run = dict(hash_run2)
    
        ## 3.b Include the controls and overlapping variants into the vcf files  
        for vcf_target in l_vcf_files_run:
    
            target_reader = vcf.Reader(filename=vcf_target)
    
            target_reader.infos['CT'] = _Info('CT',1,'String',"Controls: pool of controls that a variant appears")
            target_reader.infos['NC'] = _Info('NC',1,'Integer',"Number of Controls: number of controls that a variant appears")            
            target_reader.infos['SR'] = _Info('SR',0,'Float',"Percentage of the samples within the run with the variant")
            target_reader.infos['SR_string'] = _Info('SR_string',0,'String',"SamplesRun: pool of samples in which a variant is detected")
    
            vcf_target_annotated = os.path.splitext(vcf_target)[0]+'.annotated.vcf'
            
            l_vcf_files_run_out.append(vcf_target_annotated)
    
            target_writer = vcf.Writer(open(os.path.join(output_path,vcf_target_annotated), 'w'), target_reader)
            
            
            #hash_target  = dict(map(lambda (i,r): ((i,r.CHROM,r.POS),r), enumerate(target_reader)))
            
            hash_target = map( lambda r:(
                                        map (lambda a: ((r.CHROM,r.POS,r.REF,str(a)),r),
                                              r.ALT),
                               )[0],
                               target_reader
                              )
    
            hash_target2 = list(itertools.chain.from_iterable(hash_target))
            
            hash_target = dict(hash_target2)
            
            #for (i,chr,pos,ref,alt) in sorted(hash_target.keys()):
            for (chr,pos,ref,alt) in sorted(hash_target.keys()):
                r_target = hash_target[(chr,pos,ref,alt)]
                
                if hash_control.has_key((chr,pos,ref,alt)):
                    r_control = hash_control[(chr,pos,ref,alt)]
                    if r_control.INFO['set'] == 'Intersection':
                        r_target.INFO['CT'] = "Allcontrols_N%d" % (num_controls)
                        r_target.INFO['NC'] = num_controls
                    else:    
                        r_target.INFO['CT'] = r_control.INFO['set']
                        
                        str_controls = r_control.INFO['set']
                        
                        if str_controls == "FilteredInAll":
                            r_target.INFO['NC'] = num_controls
                        else:
                            l_ = str_controls.split('-')
                            r_target.INFO['NC'] = len(l_)
                    
                if hash_run.has_key((chr,pos,ref,alt)):
                    r_run = hash_run[(chr,pos,ref,alt)]
                    r_target.INFO['SR'] = r_run.INFO['NC']
                    # we now include the list of the samples that appears in samplesRun 
                    r_target.INFO['SR_string'] = r_run.INFO['set']

                    #if r_run.INFO['set']=='Intersection':
                    #    r_target.INFO['SR'] = True
                    
                target_writer.write_record(r_target)
                    
        sys.stdout.write('[done]\n')
        sys.stdout.flush()
    
    else:
        
        sys.stdout.write('WARNING: There are no controls included\n')
        sys.stdout.flush()
        
        ## 3.b Include only overlapping variants (controls are missing) into the vcf files  
        for vcf_target in l_vcf_files_run:
    
            target_reader = vcf.Reader(filename=vcf_target)
    
            target_reader.infos['CT'] = _Info('CT',1,'String',"Controls: pool of controls that a variant appears")
            target_reader.infos['SR'] = _Info('SR',1,'String',"Samples of the Run: The field is \'1\' if the variant appears in all the samples of the run")
    
            vcf_target_annotated = os.path.splitext(vcf_target)[0]+'.annotated.vcf'
            
            l_vcf_files_run_out.append(vcf_target_annotated)
    
            target_writer = vcf.Writer(open(os.path.join(output_path,vcf_target_annotated), 'w'), target_reader)
      
            vcf_to_iterate = walk_together(target_reader,vcf.Reader(filename=vcf_overlapping_variants_within_run)) #vcf_record_sort_key=(CHROM,POS,ALT)
    
            for [record_target,record_run] in vcf_to_iterate:
        
                if type(record_target) ==  _Record:
            
                    if type(record_run) ==  _Record:
                        record_target.INFO['SR'] = int(record_run.INFO['set']=='Intersection')
            
                    target_writer.write_record(record_target)
    
    t2 = time.time()
    sys.stdout.write("done [%0.3fs].\n\n" % (t2-t1))
    sys.stdout.flush()
    
    return l_vcf_files_run_out
        