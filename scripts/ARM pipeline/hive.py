#!/usr/bin/env python
import time
import json
import re
import pdb
import urllib,requests
import getopt, sys, csv

class api:
    def __init__(self, dna_cgi, login, pswd):
        self._dna_cgi = dna_cgi
        self._uploader_cgi = dna_cgi.replace("dna.cgi","dmUploader.cgi")
        self._login = login
        self._pswd = pswd
        self._cookie_jar = requests.cookies.RequestsCookieJar()

        self._error_submission_rgx = re.compile('^error')
        self._proc_status_rg = re.compile('"Total Progress"')
        self._seq_upload_status_rg = re.compile('"Sequence File Indexing"')
        self._uploader_success_rgx = re.compile('^Success')
        self._al_total_hits_rgx = re.compile('^\+,"total"')
 
    def __enter__(self):
        # avoid accidental auth leak
        login_dict = {'login': self._login, 'pswd': self._pswd}
        del self._login
        del self._pswd
        self.cmdr('login', login_dict)
        return self
 
    def __exit__(self, exc_type, exc_value, traceback):
        self.cmdr('logout')
        requests.get(self._dna_cgi, params=urllib.urlencode([('cmdr', 'logout')]))
 
    def cmdr(self, cmdr, param_dict=None, empty_output_is_error=False):
        if not param_dict:
            param_dict = dict()
        param_dict.update({'cmdr':cmdr,'raw':1})
        f = requests.get(self._dna_cgi, params=urllib.urlencode(param_dict),cookies=self._cookie_jar)
        f.raise_for_status()
        self._cookie_jar = f.cookies
        txt = f.text

        if not txt and empty_output_is_error:
            txt = "error: no results"
        return txt
    
    def downloadCmd(self,cmdr,filename,param_dict=None,empty_output_is_error=False):
        if not param_dict:
            param_dict = dict()
        param_dict.update({'cmdr':cmdr,'raw':1})
        r = requests.get(self._dna_cgi, params=urllib.urlencode(param_dict),cookies=self._cookie_jar)
        with open(filename, 'wb') as f:
            f.write(r.text)
            return True
        return False
            
    def propget(self, ids):
        if isinstance(ids, (list, tuple)):
            ids = ",".join(ids)
        return self.cmdr('propget', {'mode': 'json', 'ids': ids}, True)
     
    def propset(self, json_txt):
        return self.cmdr('propset2', {'mode': 'json', 'parse': json_txt})

    def waitForProcess(self,req = 0,obj = 0, status_regexp = 0):
        status = 0
        proc_status = 0
        if not status_regexp:
            status_regexp = self._proc_status_rg
        while (proc_status < 5 and status < 5):
            if obj > 0:
                res = self.cmdr('-qpRawCheck', {'raw': '1','showreqs':'0','reqObjID':obj})
            elif req > 0:
                res = self.cmdr('-qpRawCheck', {'raw': '1','showreqs':'0','req':req})
            else:
                return 0;
            for resline in res.splitlines():
                if(status_regexp.match(resline)):
                    status = int(resline.split(",")[6])
                if(self._proc_status_rg.match(resline)):
                    proc_status = int(resline.split(",")[6])
            if proc_status < 5 and status < 5:
                time.sleep(30)
        return status==5
    
    def objQry(self,qry):
        return self.cmdr('objQry',{'qry':qry})
    