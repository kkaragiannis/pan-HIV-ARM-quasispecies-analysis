#!/usr/bin/env python
import time
import json
import re
import pdb
import urllib,requests
import collections
import getopt, sys, csv
from urlparse import parse_qsl, urlparse
from email.mime.text import MIMEText
import hive

def isNumber(a):
    try:
        float(a)
        return True
    except ValueError:
        return False

class pipeline(hive.api):
    def __init__(self, dna_cgi, login, pswd):
        self._dna_cgi = dna_cgi
        self._uploader_cgi = dna_cgi.replace("dna.cgi","dmUploader.cgi")
        self._login = login
        self._pswd = pswd
        self._cookie_jar = requests.cookies.RequestsCookieJar()
        self.srrID = None
        self.denovoMinLen = 500
        self.runRow = None
        self.params = {
            "final_hexagon":{
                'ILLUMINA':{
                    'prop.svc-align-hexagon.minMatchLen.25':75,
                    'prop.svc-align-hexagon.minMatchUnit.25':1
                    },
                'LS454':{
                    'prop.svc-align-hexagon.minMatchLen.25':75,
                    'prop.svc-align-hexagon.minMatchUnit.25':0
                    },
                'PACBIO_SMRT':{
                    'prop.svc-align-hexagon.minMatchLen.25':75,
                    'prop.svc-align-hexagon.minMatchUnit.25':0
                }
            }
        }
        
        self.filename= {
            "initial_hexagon":"first_hexagon_{}",
            "initial_ref_hiveseq":"initial_ref_{}",
            "initial_ref":"initial_ref_{}.fasta",
            "denovo":"denovo_{}",
            "denovo_filtered_scaffolds":"denovo_filtered_scaffolds_{}.fasta",
            "align_denovo":"align_denovo_{}",
            "final_ref_hiveseq":"final_ref_{}",
            "final_ref":"final_ref_{}.fasta",
            "mafft":"mafft_{}",
            "final_hexagon":"final_hexagon_{}",
            "final_HH":"final_HH_{}"
            }

        self._error_submission_rgx = re.compile('^error')
        self._proc_status_rg = re.compile('"Total Progress"')
        self._seq_upload_status_rg = re.compile('"Sequence File Indexing"')
        self._uploader_success_rgx = re.compile('^Success')
        self._al_total_hits_rgx = re.compile('^\+,"total"')
 
     
    def getName(self,step):
        return self.filename[step].format(self.srrID)
 
    def isPairedEnd(self,qry_ids):
        return len(qry_ids)==2
    
    def isIllumina(self):
        return self.runRow['Platform']=="ILLUMINA"
    def is454(self):
        return self.runRow['Platform']=="LS454"
    def isPacBio(self):
        return self.runRow['Platform']=="PACBIO_SMRT"
    
    def applyPlatformSpecificParms(self,parameters,step):
        if step not in self.params:
            return parameters
        if self.runRow['Platform'] not in self.params[step]:
            return parameters
        parameters.update(self.params[step][self.runRow['Platform']])
        return parameters
        
    def getFirstHexagon(self):
        name=self.getName("initial_hexagon")
        obj_id = eval(self.objQry("alloftype('svc-align-hexagon').filter({.name=='" + name + "'})"))
        if len(obj_id):
            return obj_id[0]
        return 0
        
    def submitFirstHexagon(self,qry_ids,minMatchLen=45,maxMissQueryPercent=15,considerGoodSubAlignments=1):
        name=self.getName("initial_hexagon")
        obj_id = eval(self.objQry("alloftype('svc-align-hexagon').filter({.name=='" + name + "'})"))
        if len(obj_id):
            return obj_id[0]
        global summary
        subj_ids=eval(self.objQry("alloftype('genome',{taxonomy:'HIV subtypes'})"))
        if (not subj_ids) or not len(subj_ids):
            summary += "Failed to detect HIV representative set for initial alignment"
            return 0
        parameters={
            'svc':'dna-hexagon',
            'prop.svc-align-hexagon.alignSelector.9':'svc-align-hexagon',
            'prop.svc-align-hexagon.name.11':name,
            'prop.svc-align-hexagon.query.1':','.join(map(str,qry_ids)),
            'prop.svc-align-hexagon.submitter.14':'dna-hexagon',
            'prop.svc-align-hexagon.isbackward.7.11':1,
            'prop.svc-align-hexagon.isextendtails.7.12':0,
            'prop.svc-align-hexagon.maxExtensionGaps.7.13':0,
            'prop.svc-align-hexagon.reverseEngine.7.3':0,
            'prop.svc-align-hexagon.isoptimize.7.6':1,
            'prop.svc-align-hexagon.hashStp.7.7':'auto',
            'prop.svc-align-hexagon.hashCompileStp.7.8':'auto',
            'prop.svc-align-hexagon.looseExtenderMinimumLengthPercent.7.0':66,
            'prop.svc-align-hexagon.looseExtenderMismatchesPercent.7.4':25,
            'prop.svc-align-hexagon.alignmentEngine.7.9':1,
            'prop.svc-align-hexagon.maxHashBin.7.2':50,
            'prop.svc-align-hexagon.selfQueryPosJumpInNonPerfectAlignment.7.1':1,
            'prop.svc-align-hexagon.useRedundSim.7.5':1,
            'prop.svc-align-hexagon.computeDiagonalWidth.7.10':'auto',
            'prop.svc-align-hexagon.searchRepeatsAndTrans.5.0':0,
            'prop.svc-align-hexagon.isglobal.6.7':0,
            'prop.svc-align-hexagon.costGapNext.6.5':-4,
            'prop.svc-align-hexagon.costGapOpen.6.4':-12,
            'prop.svc-align-hexagon.costMatch.6.6':5,
            'prop.svc-align-hexagon.costMismatchNext.6.2':-6,
            'prop.svc-align-hexagon.costMismatch.6.1':-4,
            'prop.svc-align-hexagon.seed.6.3':11,
            'prop.svc-align-hexagon.allowShorterEnds.6.0':0,
            'prop.svc-align-hexagon.maxNumberQuery.24':'all',
            'prop.svc-align-hexagon.scissors.28':'hiveseq',
            'prop.svc-align-hexagon.complexityRefEntropy.3.0.0':0,
            'prop.svc-align-hexagon.complexityRefWindow.3.0.0':0,
            'prop.svc-align-hexagon.complexityEntropy.3.1.0':0,
            'prop.svc-align-hexagon.complexityWindow.3.1.0':0,
            'prop.svc-align-hexagon.keepAllMatches.33':4,
            'prop.svc-align-hexagon.maxHitsPerRead.23':200,
            'prop.svc-align-hexagon.minMatchLen.25':minMatchLen,
            'prop.svc-align-hexagon.maxMissQueryPercent.2.0':maxMissQueryPercent,
            'prop.svc-align-hexagon.considerGoodSubalignments.2.0':considerGoodSubAlignments,
            'prop.svc-align-hexagon.slice.15':'500000',
            'prop.svc-align-hexagon._type':'svc-align-hexagon',
            'prop.svc-align-hexagon.subject.16':max(subj_ids),
            'prop.svc-align-hexagon.scoreFilter.21':'None','prop.svc-align-hexagon.action.0.0':'2',
            'prop.svc-align-hexagon.doubleHash.20':'0',
            'prop.svc-align-hexagon.doubleStagePerfect.30':'1','prop.svc-align-hexagon.split.17':'query'
        }
        parameters = self.applyPlatformSpecificParms(parameters, "initial_hexagon")
        obj_id=str(self.__submitProcess(parameters).obj)
        if not int(obj_id):
            summary += "Failed to submit alignment\n"
            return 0;
        if not self.waitForProcess(obj=obj_id):
            summary += "Internal error during alignment\n"
            return 0
        
        return obj_id

    def filterReferences(self,al_id, rpkm=10000):
        al_hits=self.cmdr('alCount',{'start':0,'cnt':10,'objs':al_id})
        run_reader = csv.DictReader(al_hits.split("\n"),delimiter=",")
        hits=[]
        
        max_rpkm = 0
        for hit_row in run_reader:
            if not isNumber(hit_row["RPKM"]):
                continue
            cur_rpkm = float(hit_row["RPKM"])
            if cur_rpkm > max_rpkm:
                max_rpkm = cur_rpkm
            if( cur_rpkm > rpkm ):
                hits.append(hit_row["id"])

        if not max_rpkm:
            return hits
        if not len(hits):
            return [max_rpkm]
        return hits
    
    def getInitialRefSet(self):
        name = self.getName("initial_ref")
        obj_id = eval(self.objQry("alloftype('genome').filter({.name=='" + name + "'})"))
        if (not obj_id) or not len(obj_id):
            obj_id = eval(self.objQry("alloftype('nuc-read').filter({.name=='" + name + "'})"))
            if (not obj_id) or not len(obj_id):
                return 0
            if not self.__castObject(','.join(map(str,obj_id)),'genome'):
                summary += "Failed to cast to genome initial reference\n"
        return obj_id[0];
    
    def createRefSet(self,hits):
        name = self.getName("initial_ref_hiveseq")
        
        obj_id = eval(self.objQry("alloftype('nuc-read').filter({.name=='" + name + ".fasta'})"))
        
        if len(obj_id):
            return obj_id
        
        global summary
        subj_ids=eval(self.objQry("alloftype('genome',{taxonomy:'HIV subtypes'})"))
        if (not subj_ids) or not len(subj_ids):
            summary += "Failed to detect HIV representative set for initial alignment"
            return 0
        subj_ids=subj_ids[0]
        qrystr=""
        for hit in hits:
            qrystr+=str(subj_ids)+","+str(hit)+","+str(hit)+"\n"
        parameters = {
            'svc':'dna-hiveseq',
            'prop.svc-hiveseq.name.1':name,
            'prop.svc-hiveseq.hiveseqQry.24':qrystr,
            'prop.svc-hiveseq.AlgorithmsName.11.0':0,
            'prop.svc-hiveseq.adaptersactive.13.0':0,
            'prop.svc-hiveseq.primersactive.17.0':0,
            'prop.svc-hiveseq.QualityFilter.19':0,
            'prop.svc-hiveseq.lengthSeqFilter.16':0,
            'prop.svc-hiveseq.inputMode.25':0,
            'prop.svc-hiveseq.isFastQFile.26.0':0,
            'prop.svc-hiveseq.keepOriginalID.26.1':1,
            'prop.svc-hiveseq.submitter.29':'dna-hiveseq',
            'prop.svc-hiveseq._type':'svc-hiveseq',
            'prop.svc-hiveseq.action.0.0':'2'
        }
        obj_id=str(self.__submitProcess(parameters).obj)
        if not int(obj_id):
            summary += "Failed to submit initial reference\n"
            return 0;
        if not self.waitForProcess(obj=obj_id):
            summary += "Internal error during initial reference\n"
            return 0;
        name= self.getName("initial_ref")
        obj_id = eval(self.objQry("alloftype('nuc-read').filter({.name=='" + name + "'})"))
        if (not obj_id) or not len(obj_id):
            summary += "Failed to get hiveseq object result\n"
            return 0;
        obj_id = obj_id[0]
        if not self.__castObject(obj_id,'genome'):
            summary += "Failed to cast to genome initial reference\n"
        
        return obj_id;
    
    def runDeNovo(self,reads):
        name=self.getName("denovo")
        obj_id = eval(self.objQry("alloftype('svc-idba-hybrid').filter({.name=='" + name + "'})"))
        if len(obj_id):
            return obj_id[0]
        global summary
        reference=eval(self.objQry("alloftype('genome',{taxonomy:'HXB2'})"))
        if (not reference) or not len(reference):
            summary += "Failed to detect HXB2 reference for denovo hybrid"
            return 0
        reference = reference[0]
        
        if not hasattr(self,'denovoscript'):
            self.denovoscript = eval(self.objQry("alloftype('algorithm-script',{target_page:'dna-pipeline-idba-hybrid'})[0]"))
        
        ispaired = self.isPairedEnd(reads)
        parameters={
            'svc':'generic-launcher',
            'prop.svc-idba-hybrid.name.1':name,
            'prop.svc-idba-hybrid.resubmitMode.24':0,
            'prop.svc-idba-hybrid.reads.26.0':reads[0],
            'prop.svc-idba-hybrid.reads.26.1':(reads[1] if ispaired else reads[0]),
            'prop.svc-idba-hybrid.mink.27':20,
            'prop.svc-idba-hybrid.maxk.28':124,
            'prop.svc-idba-hybrid.step.29':10,
            'prop.svc-idba-hybrid.prefix.30':3,
            'prop.svc-idba-hybrid.min_count.31':2,
            'prop.svc-idba-hybrid.min_support.32':1,
            'prop.svc-idba-hybrid.seed_kmer.33':11,
            'prop.svc-idba-hybrid.min_contig.34':200,
            'prop.svc-idba-hybrid.similar.35':0.80,
            'prop.svc-idba-hybrid.max_mismatch.36':3,
            'prop.svc-idba-hybrid.min_pairs.37':3,
            'prop.svc-idba-hybrid.submitter.42':'dna-pipeline-idba-hybrid',
            'prop.svc-idba-hybrid.max_gap.14':500,
            'prop.svc-idba-hybrid.min_region.15':'3000',
            'prop.svc-idba-hybrid._type':'svc-idba-hybrid',
            'prop.svc-idba-hybrid.reference.19':reference,
            'prop.svc-idba-hybrid.algo.6':self.denovoscript,
            'prop.svc-idba-hybrid.service.21':'generic-launcher',
            'prop.svc-idba-hybrid.nrepeat.22.0':1,
            'prop.svc-idba-hybrid.svcTitle.0.0':'Generic External Application Launcher',
            'prop.svc-idba-hybrid.action.0.0':2,
            'prop.svc-idba-hybrid.folder.0.3':'Inbox',
            'prop.svc-idba-hybrid.isPostponed.23.3':0,
            'prop.svc-idba-hybrid.reqPriority.23.7':0
        }
        parameters = self.applyPlatformSpecificParms(parameters, "denovo")
        obj_id=str(self.__submitProcess(parameters).obj)
        if not int(obj_id):
            summary += "Failed to submit DeNovo\n"
            return 0;
        if not self.waitForProcess(obj=obj_id):
            summary += "Internal error during DeNovo\n"
            return 0;
        return obj_id
    
    def uploaded_denovo_results(self,filename):
        return eval(self.objQry("alloftype('nuc-read').filter({.name=~/" + filename + "/})"))
        
    def getDeNovo_results(self,denovoID,filename):
        params={
            'ids':denovoID,
            'filename':'contig.fa'
        }
        tmp_filename=filename+"_tmp"
        self.downloadCmd("objFile",tmp_filename,params)
                
        if (not os.path.isfile(tmp_filename)) or (os.stat(tmp_filename).st_size == 0 ) :
            if (os.path.isfile(tmp_filename)):
                os.remove(tmp_filename)
            return False
    
        seqs=[]
        for record in SeqIO.parse(tmp_filename, "fasta"):
            if len(record.seq) > self.denovoMinLen:
                seqs.append(record)
        if not len(seqs):
            if os.path.isfile(tmp_filename):
                os.remove(tmp_filename)
            return False
        SeqIO.write(seqs,filename, "fasta")
        if os.path.isfile(tmp_filename):
            os.remove(tmp_filename)
        
        return True
            
    def uploadDeNovo_results(self,filename):
        global summary
        upl_id = self.uploadFile(filename)
        if not self.waitForProcess(obj=upl_id) :
            summary += "Something went wrong during processing of uploaded file\n"
            return 0
        ref_id = eval(self.objQry("alloftype('nuc-read').filter({.name=~/" + filename + "/})"))
        if not len(ref_id):
            summary += "Something went wrong during denovo filtered scaffolds upload\n"
            return 0
        
        os.remove(filename)
        return ref_id
    
    def getAlignDenovoResults(self):
        name=self.getName("align_denovo")
        obj_id = eval(self.objQry("alloftype('svc-align-hexagon').filter({.name=='" + name + "'})"))
        if len(obj_id):
            return obj_id[0]
        return 0
    
    def alignDeNovoResult(self, scaffoldID):
        name=self.getName("align_denovo")
        
        obj_id = eval(self.objQry("alloftype('svc-align-hexagon').filter({.name=='" + name + "'})"))
        if len(obj_id):
            return obj_id[0]
        
        refID=eval(self.objQry("alloftype('genome',{taxonomy:'HXB2'})"))
        
        parameters={
            'svc':'dna-hexagon',
            'prop.svc-align-hexagon.alignSelector.9':'svc-align-hexagon',
            'prop.svc-align-hexagon.name.11':name,
            'prop.svc-align-hexagon.query.1':','.join(map(str,scaffoldID)),
            'prop.svc-align-hexagon.submitter.14':'dna-hexagon',
            'prop.svc-align-hexagon.isbackward.7.11':1,
            'prop.svc-align-hexagon.isextendtails.7.12':0,
            'prop.svc-align-hexagon.maxExtensionGaps.7.13':0,
            'prop.svc-align-hexagon.reverseEngine.7.3':0,
            'prop.svc-align-hexagon.isoptimize.7.6':1,
            'prop.svc-align-hexagon.hashStp.7.7':'auto',
            'prop.svc-align-hexagon.hashCompileStp.7.8':'auto',
            'prop.svc-align-hexagon.looseExtenderMinimumLengthPercent.7.0':66,
            'prop.svc-align-hexagon.looseExtenderMismatchesPercent.7.4':25,
            'prop.svc-align-hexagon.alignmentEngine.7.9':1,
            'prop.svc-align-hexagon.maxHashBin.7.2':50,
            'prop.svc-align-hexagon.selfQueryPosJumpInNonPerfectAlignment.7.1':1,
            'prop.svc-align-hexagon.useRedundSim.7.5':1,
            'prop.svc-align-hexagon.computeDiagonalWidth.7.10':60,
            'prop.svc-align-hexagon.searchRepeatsAndTrans.5.0':0,
            'prop.svc-align-hexagon.isglobal.6.7':0,
            'prop.svc-align-hexagon.costGapNext.6.5':-4,
            'prop.svc-align-hexagon.costGapOpen.6.4':-12,
            'prop.svc-align-hexagon.costMatch.6.6':5,
            'prop.svc-align-hexagon.costMismatchNext.6.2':-6,
            'prop.svc-align-hexagon.costMismatch.6.1':-4,
            'prop.svc-align-hexagon.seed.6.3':11,
            'prop.svc-align-hexagon.allowShorterEnds.6.0':0,
            'prop.svc-align-hexagon.maxNumberQuery.24':'all',
            'prop.svc-align-hexagon.scissors.28':'hiveseq',
            'prop.svc-align-hexagon.complexityRefEntropy.3.0.0':0,
            'prop.svc-align-hexagon.complexityRefWindow.3.0.0':0,
            'prop.svc-align-hexagon.complexityEntropy.3.1.0':0,
            'prop.svc-align-hexagon.complexityWindow.3.1.0':0,
            'prop.svc-align-hexagon.keepAllMatches.33':4,
            'prop.svc-align-hexagon.maxHitsPerRead.23':200,
            'prop.svc-align-hexagon.minMatchLen.25':75,
            'prop.svc-align-hexagon.maxMissQueryPercent.2.0':45,
            'prop.svc-align-hexagon.considerGoodSubalignments.2.0':1,
            'prop.svc-align-hexagon.slice.15':'500000',
            'prop.svc-align-hexagon._type':'svc-align-hexagon',
            'prop.svc-align-hexagon.subject.16':max(refID),
            'prop.svc-align-hexagon.scoreFilter.21':'None','prop.svc-align-hexagon.action.0.0':'2',
            'prop.svc-align-hexagon.doubleHash.20':'0',
            'prop.svc-align-hexagon.doubleStagePerfect.30':'1','prop.svc-align-hexagon.split.17':'query'}
        
        parameters = self.applyPlatformSpecificParms(parameters, "align_denovo")
        obj_id=str(self.__submitProcess(parameters).obj)
        global summary
        if not int(obj_id):
            summary += "Failed to submit alignment\n"
            return 0;
        if not self.waitForProcess(obj=obj_id):
            summary += "Internal error during alignment\n"
            return 0
        
#         if not self.waitForProcess(obj=obj_id,status_regexp=self._seq_upload_status_rg) :
#             summary += "Something went wrong during first hexagon\n"
#             return
        return obj_id
                
    def getValidatedDeNovoReference(self,alid):
        params = {
            'start':0,
            'mySubID':1,
            'objs':alid,
            'found':1
        }
        response = self.cmdr("alMatch", params)
        reader = csv.DictReader(response.split("\n"),delimiter=",")
        
        directionalCoverage=[{'ids':[],'totalLength':0,'ranges':[]},
                             {'ids':[],'totalLength':0,'ranges':[]}]
        for crow in reader:
            ref_start = int(crow['Reference Start'])
            ref_end = int(crow['Reference End'])
            length = ref_end - ref_start
            cdir = int(crow['Direction'])
            if length < self.denovoMinLen:
                continue;
            if cdir==-1 :
                cdir_ind = 0
            else:
                cdir_ind = 1
            ldir = directionalCoverage[cdir_ind]
            t_start = ref_start
            t_end = ref_end
            for crng in ldir['ranges']:
                if crng['start'] == crng['end'] or crng['end'] < t_start or t_end < crng['start']:
                    continue
                if crng['start'] > t_start and crng['end'] < t_end:
                    ldir.totlalength -= crng['end'] - crng['start']
                    crng['end'] = crng['start']
                elif crng['start'] < t_start and crng['end'] > t_end:
                    t_start = t_end
                elif crng['start'] < t_start and crng['end'] < t_end:
                    t_start = crng['end']
                else:
                    t_end = crng['start']
            if t_end > t_start:
                ldir['totalLength'] += t_end - t_start
                ldir['ranges'].append({'start':t_start,'end':t_end})
                ldir['ids'].append(crow['Read #'])
        max_i = 1
        if directionalCoverage[0]['totalLength'] > directionalCoverage[1]:
            max_i = 0
        
        return directionalCoverage[max_i]['ids']
            
    def getCreatedReferenceSequence(self):
        global summary
        name=self.getName("final_ref")
        obj_id  = eval(self.objQry("alloftype('nuc-read').filter({.name=~/" + name + "/})"))
        if obj_id and len(obj_id):
            obj_id = obj_id[0]
            if not self.__castObject(obj_id,'genome'):
                summary += "Failed to cast to genome final reference set\n"
                return 0
            return obj_id
        obj_id = eval(self.objQry("alloftype('genome').filter({.name=~/" + name + "/})"))
        if obj_id and len(obj_id):
            return obj_id[0]
        return 0
        

    def createReferenceSequence(self,denovoRefID,denovoIDlist,refID):
        qryStr = ""
        for crid in denovoIDlist:
            qryStr += str(denovoRefID)+","+str(crid)+","+str(crid)+"\n"
        qryStr+=str(refID)+",0,100000\n"
        
        name=self.getName("final_ref_hiveseq")
        
        parameters = {
            'svc':'dna-hiveseq',
            'prop.svc-hiveseq.name.1':name,
            'prop.svc-hiveseq.hiveseqQry.24':qryStr,
            'prop.svc-hiveseq.AlgorithmsName.11.0':0,
            'prop.svc-hiveseq.adaptersactive.13.0':0,
            'prop.svc-hiveseq.primersactive.17.0':0,
            'prop.svc-hiveseq.QualityFilter.19':0,
            'prop.svc-hiveseq.lengthSeqFilter.16':0,
            'prop.svc-hiveseq.inputMode.25':0,
            'prop.svc-hiveseq.isFastQFile.26.0':0,
            'prop.svc-hiveseq.keepOriginalID.26.1':1,
            'prop.svc-hiveseq.submitter.29':'dna-hiveseq',
            'prop.svc-hiveseq._type':'svc-hiveseq',
            'prop.svc-hiveseq.action.0.0':'2'
        }
        obj_id=str(self.__submitProcess(parameters).obj)
        global summary
        if not int(obj_id):
            summary += "Failed to submit final reference creation\n"
            return 0
        if not self.waitForProcess(obj=obj_id):
            summary += "Internal error during final reference creation\n"
            return 0
        name=self.getName("final_ref")
        obj_id = eval(self.objQry("alloftype('nuc-read').filter({.name=='" + name + "'})"))
        if not len(obj_id):
            summary += "Failed to get final reference creation\n"
            return 0
        obj_id = obj_id[0]
        if not self.__castObject(obj_id,'genome'):
            summary += "Failed to cast to genome final reference set\n"
            return 0
        
        return obj_id;
    
    def getMAFFT(self):
        name=self.getName("mafft")
        
        obj_id  = eval(self.objQry("alloftype('svc-align-mafft').filter({.name=~/" + name + "/})"))
        if obj_id and len(obj_id):
            return obj_id[0]
        return 0
    
    def runMAFFT(self,refID):
        name=self.getName("mafft")
        
        obj_id  = eval(self.objQry("alloftype('svc-align-mafft').filter({.name=~/" + name + "/})"))
        if obj_id and len(obj_id):
            return obj_id[0]
        dim = eval(self.objQry("((obj)"+str(refID)+")['rec-count']"))
        if dim <= 1:
            return 0
        parameters = {
            "svc":"dna-alignx",
            "prop.svc-align-mafft.name.19":name,
            "prop.svc-align-mafft.alignSelector.18":"svc-align-mafft",
            "prop.svc-align-mafft.submitter.20":"dna-hexagon&cmdMode=mafft",
            "prop.svc-align-mafft.scissors.12":"hiveseq",
            "prop.svc-align-mafft._type":"svc-align-mafft",
            "prop.svc-align-mafft.random_seed.11":0,
            "prop.svc-align-mafft.subject.15":refID,
            "prop.svc-align-mafft.slice.14.1":500000,
            "prop.svc-align-mafft.action.17.0":2,
            "prop.svc-align-mafft.folder.17.2":"Inbox",
            "prop.svc-align-mafft.reqPriority.17.7":0,
            "prop.svc-align-mafft.split.13":"query"
        }
        obj_id=str(self.__submitProcess(parameters).obj)
        global summary
        if not int(obj_id):
            summary += "Failed to submit MAFFT\n"
            return 0;
        if not self.waitForProcess(obj=obj_id):
            summary += "Internal error during MAFFT\n"
            return 0;
        return obj_id

    def getFinalHexagon(self):
        name = self.getName("final_hexagon")
        
        obj_id  = eval(self.objQry("alloftype('svc-align-hexagon').filter({.name=~/" + name + "/})"))
        if obj_id and len(obj_id):
            return obj_id[0]
        return 0
    
    def submitFinalHexagon(self,qry_ids,ref_id):
        name = self.getName("final_hexagon")
        
        obj_id  = eval(self.objQry("alloftype('svc-align-hexagon').filter({.name=~/" + name + "/})"))
        if obj_id and len(obj_id):
            return obj_id[0]
        
        parameters={
            'svc':'dna-hexagon',
            'prop.svc-align-hexagon.alignSelector.9':'svc-align-hexagon',
            'prop.svc-align-hexagon.name.11':name,
            'prop.svc-align-hexagon.query.1':','.join(map(str,qry_ids)),
            'prop.svc-align-hexagon.submitter.14':'dna-hexagon',
            'prop.svc-align-hexagon.isbackward.7.11':1,
            'prop.svc-align-hexagon.isextendtails.7.12':0,
            'prop.svc-align-hexagon.maxExtensionGaps.7.13':0,
            'prop.svc-align-hexagon.reverseEngine.7.3':0,
            'prop.svc-align-hexagon.isoptimize.7.6':1,
            'prop.svc-align-hexagon.hashStp.7.7':'auto',
            'prop.svc-align-hexagon.hashCompileStp.7.8':'auto',
            'prop.svc-align-hexagon.looseExtenderMinimumLengthPercent.7.0':66,
            'prop.svc-align-hexagon.looseExtenderMismatchesPercent.7.4':25,
            'prop.svc-align-hexagon.alignmentEngine.7.9':1,
            'prop.svc-align-hexagon.maxHashBin.7.2':50,
            'prop.svc-align-hexagon.selfQueryPosJumpInNonPerfectAlignment.7.1':1,
            'prop.svc-align-hexagon.useRedundSim.7.5':1,
            'prop.svc-align-hexagon.computeDiagonalWidth.7.10':60,
            'prop.svc-align-hexagon.searchRepeatsAndTrans.5.0':0,
            'prop.svc-align-hexagon.isglobal.6.7':0,
            'prop.svc-align-hexagon.costGapNext.6.5':-4,
            'prop.svc-align-hexagon.costGapOpen.6.4':-12,
            'prop.svc-align-hexagon.costMatch.6.6':5,
            'prop.svc-align-hexagon.costMismatchNext.6.2':-6,
            'prop.svc-align-hexagon.costMismatch.6.1':-4,
            'prop.svc-align-hexagon.seed.6.3':11,
            'prop.svc-align-hexagon.allowShorterEnds.6.0':0,
            'prop.svc-align-hexagon.maxNumberQuery.24':'all',
            'prop.svc-align-hexagon.scissors.28':'hiveseq',
            'prop.svc-align-hexagon.complexityRefEntropy.3.0.0':0,
            'prop.svc-align-hexagon.complexityRefWindow.3.0.0':0,
            'prop.svc-align-hexagon.complexityEntropy.3.1.0':0,
            'prop.svc-align-hexagon.complexityWindow.3.1.0':0,
            'prop.svc-align-hexagon.keepAllMatches.33':4,
            'prop.svc-align-hexagon.maxHitsPerRead.23':200,
            'prop.svc-align-hexagon.minMatchLen.25':65,
            'prop.svc-align-hexagon.minMatchUnit.25': 1,
            'prop.svc-align-hexagon.trimLowScoreEndsMaxMismatches.80.0': 10,
            'prop.svc-align-hexagon.trimLowScoreEndsWindow.80.0': 20,
            'prop.svc-align-hexagon.maxMissQueryPercent.2.0':45,
            'prop.svc-align-hexagon.considerGoodSubalignments.2.0':1,
            'prop.svc-align-hexagon.slice.15':'500000',
            'prop.svc-align-hexagon._type':'svc-align-hexagon',
            'prop.svc-align-hexagon.subject.16':ref_id,
            'prop.svc-align-hexagon.scoreFilter.21':'None',
            'prop.svc-align-hexagon.action.0.0':'2',
            'prop.svc-align-hexagon.doubleHash.20':'0',
            'prop.svc-align-hexagon.doubleStagePerfect.30':'1',
            'prop.svc-align-hexagon.split.17':'query'}
        parameters = self.applyPlatformSpecificParms(parameters, "final_hexagon")
        obj_id=str(self.__submitProcess(parameters).obj)
        global summary
        if not int(obj_id):
            summary += "Failed to submit alignment\n"
            return 0;
        if not self.waitForProcess(obj=obj_id):
            summary += "Internal error during alignment\n"
            return 0
        
        return obj_id

    def getHexaHedron(self):
        name=self.getName("final_HH")
        obj_id  = eval(self.objQry("alloftype('svc-popul').filter({.name=~/" + name + "/})"))
        if obj_id and len(obj_id):
            return obj_id[0]
        return 0        

    def submitHexaHedron(self,al_id,mafft_ids,ispaired,mut_threshold=5):
        name=self.getName("final_HH")
        obj_id  = eval(self.objQry("alloftype('svc-popul').filter({.name=~/" + name + "/})"))
        if obj_id and len(obj_id):
            return obj_id[0]
        parameters = {
            'svc':'dna-popul',
            'prop.svc-popul.sortClones.5':'0',
            'prop.svc-popul.pairedEnd.7':ispaired,
            'prop.svc-popul.fuzzy_threshold.8':'20',
            'prop.svc-popul.name.9':name,
            'prop.svc-popul.treeSize.1.5':'500',
            'prop.svc-popul.TreeFrameSize.1.0':'30',
            'prop.svc-popul.TreeStep.1.1':'30',
            'prop.svc-popul.TreeReadOverlap.1.2':'1',
            'prop.svc-popul.InPopulationFrameTrees.1.3':'0',
            'prop.svc-popul.generateTrees.1.4':'0',
            'prop.svc-popul.submitter.10':'dna-popul',
            'prop.svc-popul.parent_proc_ids.12':al_id,
            'prop.svc-popul._type':'svc-popul',
            'prop.svc-popul.action.0.0':'2',
            'prop.svc-popul.bifThreshold.15':mut_threshold
        }
        if mafft_ids and len(mafft_ids):
            parameters['prop.svc-popul.mutualAligmentID.3.2'] = mafft_ids
            
        obj_id=str(self.__submitProcess(parameters).obj)
        global summary
        if not int(obj_id):
            summary += "Failed to submit HexaHedron\n"
            return 0;
        if not self.waitForProcess(obj=obj_id):
            summary += "Internal error during HexaHedron\n"
            return 0;
        return obj_id
    
    def initRun(self,row):
        self.srrID = row['Run']
        self.runRow = row
    
    def singleRun(self):
        global summary
        
        srrID = self.srrID

        summary = "Run: "+srrID+"\n"
        if self.getHexaHedron():
            return 1
        
        qry_ids = eval(self.objQry("alloftype('nuc-read').filter({.name=~/" + srrID + "/})"))
        ispaired = self.isPairedEnd(qry_ids)
        
        if not qry_ids:
            summary += "Cannot find reads"
            return 0
        final_alid = self.getFinalHexagon()
        mafft_id=''
        if (not final_alid):
            mafft_id = self.getMAFFT()
            if (not mafft_id):
                final_refID = self.getCreatedReferenceSequence()
                if (not final_refID):
                    denovo_alid = self.getAlignDenovoResults()
                    refSet_ID=""
                    denovo_filtered_reflist=[]
                    denovo_scaff_id = ""
                    filename=self.getName("denovo_filtered_scaffolds")
                    if (not denovo_alid) :
                        refSet_ID = self.getInitialRefSet()
                        if (not refSet_ID):
                            al_id=self.submitFirstHexagon(qry_ids)
                            if not al_id:
                                return 0
                            
                            filtered_hits=self.filterReferences(al_id)
                            if not len(filtered_hits):
                                summary += "No hits generated from first hexagon"
                                return 1
        
                            refSet_ID = self.createRefSet(filtered_hits)
                            if not refSet_ID:
                                return 0
                            
                        denovo_ID = self.runDeNovo(qry_ids)
                        if not denovo_ID:
                            return 0
                        denovo_scaff_id = self.uploaded_denovo_results(filename)
                        if (not denovo_scaff_id) or not len(denovo_scaff_id):
                            if (not self.getDeNovo_results(denovo_ID, filename)) or (not os.path.isfile(filename)) or (os.stat(filename).st_size == 0 ):
                                summary += "Failed to generate denovo assembly scaffolds file"
                            else:
                                denovo_scaff_id = self.uploadDeNovo_results(filename)
                                if not denovo_scaff_id:
                                    return 0
                        if len(denovo_scaff_id):
                            denovo_alid = self.alignDeNovoResult(denovo_scaff_id)
                            if not denovo_alid:
                                return 0
                    else:
                        refSet_ID = self.getInitialRefSet()
                        denovo_scaff_id = self.uploaded_denovo_results(filename)
                    if os.path.isfile(filename):
                        os.remove(filename)
                    if denovo_alid:
                        denovo_filtered_reflist = self.getValidatedDeNovoReference(denovo_alid)
                    
                    if not refSet_ID:
                        summary += "No initial RefID\n"
                        return 0                   
                    if not len(denovo_scaff_id):
                        summary += "No valid scaffolds generated\n"
                    else:
                        denovo_scaff_id = denovo_scaff_id[0]
                    if not len(denovo_filtered_reflist):
                        summary += "All scaffolds were filtered out\n"
                        
                    final_refID = self.createReferenceSequence(denovo_scaff_id, denovo_filtered_reflist, refSet_ID)
                    if not final_refID:
                        return 0
                mafft_id = self.runMAFFT(final_refID)
            
            final_alid = self.submitFinalHexagon(qry_ids, final_refID)
        if not final_alid:
            return 0
        
        hxh_id=self.submitHexaHedron(final_alid,mafft_id,ispaired)
        if not hxh_id:
            summary += "Failed to submit main HexaHedron job"
            return 0
        return 1
    
    def getUnaligned(self,isInitial):
        al_id = None
        if isInitial:
            al_id = self.getFirstHexagon()
        else:
            al_id = self.getFinalHexagon()
        if not al_id:
            return 0
        
        al_hits=self.cmdr('alCount',{'start':0,'cnt':10,'objs':al_id})
        run_reader = csv.DictReader(al_hits.split("\n"),delimiter=",")
        
        for hit_row in run_reader:
            if hit_row["Reference"]=="Unaligned":
                return hit_row["Hits Unique"]
        return 0
    
#     def getRunStats(self):
#         qry_ids = eval(self.objQry("alloftype('nuc-read').filter({.name=~/^"+self.srrID+"/})"))
#         ispaired = self.isPairedEnd(qry_ids)
#         
#         read_cnt = eval(self.objQry("("+str(qry_ids)+" as objlist).map({this['rec-count']}).sum()"))
# 
#         init_unalign_count = self.getUnaligned(True)
#         
#         init_ref_count = 0
#         init_ref_ids = eval(self.objQry("alloftype('nuc-read').filter({.name=='"+self.getName("initial_ref")+"'})"))
#         if init_ref_ids and len(init_ref_ids):
#             init_ref_count = eval(self.objQry("("+str(init_ref_ids)+" as objlist)[0]['rec-count']"))
#         
#         contig_count = 0
#         contig_ids = eval(self.objQry("alloftype('nuc-read').filter({.name=='"+self.getName("denovo_filtered_scaffolds")+"'})"))
#         if contig_ids and len(contig_ids):
#             contig_count = eval(self.objQry("("+str(contig_ids)+" as objlist)[0]['rec-count']"))
# 
#         final_ref_count = 0
#         final_ref_ids = eval(self.objQry("alloftype('genome').filter({.name=='"+self.getName("final_ref")+"'})"))
#         if final_ref_ids and len(final_ref_ids):
#             final_ref_count = eval(self.objQry("("+str(final_ref_ids)+" as objlist)[0]['rec-count']"))
#                     
#         final_unalign_count = self.getUnaligned(False)
#         
#         return [ispaired,read_cnt,init_unalign_count,init_ref_count,
#                 contig_count,final_ref_count,final_unalign_count]
#         
#         
# def getRunStats(hexahedronWF):
#     qry_ids = eval(hexahedronWF.objQry("alloftype('nuc-read').filter({.name=~/^"+hexahedronWF.srrID+"/})"))
#     ispaired = hexahedronWF.isPairedEnd(qry_ids)
#     
#     read_cnt = eval(hexahedronWF.objQry("("+str(qry_ids)+" as objlist).map({this['rec-count']}).sum()"))
# 
#     init_unalign_count = hexahedronWF.getUnaligned(True)
#     
#     init_ref_count = 0
#     init_ref_ids = eval(hexahedronWF.objQry("alloftype('genome').filter({.name=='"+hexahedronWF.getName("initial_ref")+"'})"))
#     if init_ref_ids and len(init_ref_ids):
#         init_ref_count = eval(hexahedronWF.objQry("("+str(init_ref_ids)+" as objlist)[0]['rec-count']"))
#     
#     contig_count = 0
#     contig_ids = eval(hexahedronWF.objQry("alloftype('nuc-read').filter({.name=='"+hexahedronWF.getName("denovo_filtered_scaffolds")+"'})"))
#     if contig_ids and len(contig_ids):
#         contig_count = eval(hexahedronWF.objQry("("+str(contig_ids)+" as objlist)[0]['rec-count']"))
# 
#     final_ref_count = 0
#     final_ref_ids = eval(hexahedronWF.objQry("alloftype('genome').filter({.name=='"+hexahedronWF.getName("final_ref")+"'})"))
#     if final_ref_ids and len(final_ref_ids):
#         final_ref_count = eval(hexahedronWF.objQry("("+str(final_ref_ids)+" as objlist)[0]['rec-count']"))
#                 
#     final_unalign_count = hexahedronWF.getUnaligned(False)
#     
#     return [ispaired,read_cnt,init_unalign_count,init_ref_count,
#             contig_count,final_ref_count,final_unalign_count]