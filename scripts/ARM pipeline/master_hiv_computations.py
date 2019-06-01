#!/usr/bin/env python
import time
import json, csv
import re
import urllib,requests
import getpass
import collections
import getopt, sys
import smtplib
import os
import shutil
import pdb
from Bio import SeqIO
from email.MIMEText import MIMEText
from multiprocessing import Pool
import hivehexahedronworkflow as hhwf


hive_deployment = 'http://hive.biochemisty.gwu.edu/dna.cgi'

def runSamples(hive,inputfilepath,progressfilepath,chnk,chnkcnt):
    prgr_header = [
        "#","srrID","Status","Paired End","Total Reads","Unaligned (qual)",
        "Reference #","Assembled contigs","Final Reference #", "Unaligned (final)"
        ]
    
    statusfilepath=progressfilepath+"_tmp"
    scannedsrrID = {}
    if os.path.isfile(progressfilepath):
        with open(progressfilepath, 'rb') as prgfile:
            prg_reader = csv.DictReader(prgfile)
            for run_row in prg_reader:
                scannedsrrID[run_row['srrID']] = True
        
    global summary
    summary=""
    with open(progressfilepath, 'a+',0) as prgfile, open(statusfilepath,'w+',0) as statusfile:
        prg_writer = csv.writer(prgfile)
        
        if ( not os.path.isfile(progressfilepath) ) or os.stat(progressfilepath).st_size == 0 :
            prg_writer.writerow(prgr_header)
        
        cnt=len(scannedsrrID.keys())
        srrID=0
        error_msg=True
        irow = -1
        with open(inputfilepath) as inputfile:
            in_reader = csv.DictReader(inputfile)
            summary =""

            for run_row in in_reader:
                srrID = run_row["Run"]
                irow += 1
                if irow % chnkcnt != chnk:
                    continue
                
                if srrID in scannedsrrID:
                    statusfile.write(srrID+": processed already\n")
                    continue
                
                hive.initRun(run_row)
                is_success = hive.singleRun() 
                if not is_success:
                    prgrRecord = [cnt,srrID,'failed']
                    statusfile.write("Failed\n")
                    error_msg = True
#                     break
                else:
                    prgrRecord = [cnt,srrID,'completed']
                    prgrRecord.extend(hhwf.getRunStats(hive))
                    statusfile.write("Completed\n")
                statusfile.write(summary+"\n")
                cnt += 1
                
                
                prg_writer.writerow(prgrRecord)
    return 1
#     email_body=""
#     if error_msg:
#         email_body +='Error in Run <b>'+str(srrID)+'</b><br>Summary:<br>'+summary
#     else:
#         email_body +='Success'
    
#     msg = MIMEText(email_body,'html')
#     dst_mail = '*****'
#     msg['Subject'] = 'Your HexaHedron pipeline update'
#     msg['From'] = '****'
#     msg['To'] = dst_mail
#      
#     s = smtplib.SMTP('localhost')
#     s.sendmail('****', [dst_mail], msg.as_string())
#     s.quit()

        
def executefunc(args):
    runSamples( *args )

def main() :
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'i:o:s:u:p:h:', ['input=', 'output=', 'status=', 'user=', 'password=', 'help'])
    except getopt.GetoptError as err:
        print str(err)
        usage()
        sys.exit(2)
    
    for opt, arg in opts:
        if opt in ('-h', '--help'):
            usage()
            sys.exit(2)
        elif opt in ('-i', '--input'):
            inputfilepath = arg
        elif opt in ('-p', '--progress'):
            progressfilepath = arg
        elif opt in ('-u', '--user'):
            user = arg
            pswd = getpass.getpass('Password:')
        else:
            usage()
            sys.exit(2)

    global hive_deployment
    with hhwf.pipeline(hive_deployment, user, pswd) as hive:
        chnkcnt = 32
        chnk_progressfilepath = [(progressfilepath+"chnk_{}").format(i) for i in range(chnkcnt)]
        chnk_args = [[hive,inputfilepath,chnk_progressfilepath[i],i,chnkcnt] for i in range(chnkcnt)]
        pool = Pool(chnkcnt)
        pool.map(executefunc,chnk_args)

#         chnkcnt = 1
#         chnk_progressfilepath = [(progressfilepath+"chnk_{}").format(i) for i in range(chnkcnt)]
#         chnk_args = [[hive,inputfilepath,chnk_progressfilepath[i],i,chnkcnt] for i in range(chnkcnt)]
#         executefunc(chnk_args[0])


if __name__ == '__main__':
    main()

