#!/usr/bin/env python
import time
import json
import re, os, csv, random
import pdb
import urllib,requests
import collections
import getopt, sys
from Bio import SeqIO
from gql import gql
from sierrapy import SierraClient

results_destinations = "E:\\app\\HexaHedron_scripts\\pan_HIV_analysis\\results\\pipeline_sequences\\"

rs_query = """query example($sequences: [UnalignedSequenceInput]!) {
    viewer {
        currentVersion {
          text
          publishDate
        }
        sequenceAnalysis(sequences: $sequences) {
          validationResults {
            level
            message
          }
          inputSequence {
            header
          },
          drugResistance {
            drugScores {
              drugClass {
                name
                fullName
              }
              drug {
                displayAbbr
              }
              SIR
              score
              level
              text
              partialScores {
                mutations {
                  text
                }
                score
              }
            }
          }
        }
      }
    }"""


def runHIVDBanalysis(flnm):
    
    if not os.path.exists(flnm):
        return ""
    
    sequences = []
    for rec in SeqIO.parse(flnm,'fasta'):
        sequences.extend([
            {
                'header':rec.id,
                'sequence':str(rec.seq)
            }]
        )
        
    global rs_query
    
    if not sequences or not len(sequences):
        return ""

    client = SierraClient('https://hivdb.stanford.edu/graphql')
    
    data = client.execute( gql(rs_query), variable_values={"sequences":sequences})
    
    return data
    
def printHIVDBresults(results, srrID, assembly, rec_cnt, out_fh):
    
    if not results:
        rec_cnt += 1
        row = {}
        row["srrID"] = srrID
        row["assembly"] = assembly
        row["#"] = rec_cnt
        row["sequence"] = "missing"
        out_fh.writerow(row)
        return rec_cnt
    
    sequenceAnalysis = results["viewer"]["sequenceAnalysis"]
    cnt_rpts = 0
    rows = []
    for seqs in sequenceAnalysis:
        drugResArr = seqs["drugResistance"]
        for drugRes in drugResArr:
            for drugScore in drugRes["drugScores"]:
                row = {}
                row["srrID"] = srrID
                row["assembly"] = assembly
                row["sequence"] = seqs["inputSequence"]["header"]
                row["class"] = drugScore["drugClass"]["name"]
                row["drug"] = drugScore["drug"]["displayAbbr"]
                row["score"] = drugScore["score"]
                row["level"] = drugScore["level"]
                mutsjson = {}
                mutsjson["mutations_arr"] = []
                mutsjson["mutations_arr"]
                my_mut = []
                for mutRes in drugScore["partialScores"]:
                    my_mut = []
                    for muts in mutRes["mutations"]:
                        my_mut.append(muts["text"])
                    score = mutRes["score"]
                    
                    mutsjson["mutations_arr"].append({'mutations':my_mut,'score':score})
                
                if len(mutsjson["mutations_arr"]) == 0 :
                    continue
                row["mutations"] = json.dumps(mutsjson["mutations_arr"])
                rec_cnt += 1
                cnt_rpts += 1
                row["#"] = rec_cnt
                #this way everything from one srrID-assembly is guaranteed to be there and we can resume if something was stopped
                rows.append(row)

    if cnt_rpts == 0:
        rec_cnt += 1
        row = {}
        row["srrID"] = srrID
        row["assembly"] = assembly
        row["#"] = rec_cnt
        row["sequence"] = "-"
        out_fh.writerow(row)
    else:
        #this way everything from one srrID-assembly is guaranteed to be there and we can resume if something was stopped
        for row in rows:
            out_fh.writerow(row)     
    return rec_cnt
   
def usage():
   filename = os.path.basename(__file__)
   print """
   usage: {file} [OPTIONS]
   
   options:
      -h, --help, -?                Print this usage info
      -i, --input                   path to the hxh_progress.csv file that has the path of all fasta files
      -o, --output                  directory output where all the results are in (default current directory)
   """.format(file=filename)


def get_path(input_path):
    if pd.isna(input_path):
        return ""
    filename = input_path.split("/")[-1]
    global results_destinations
    return results_destinations+filename


def main() :
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'i:o:h', ['input=', 'output=', 'help'])
    except getopt.GetoptError as err:
        print str(err)
        usage()
        sys.exit(2)
    
    output = ""
    input = ""
    
    header = ["#","srrID","assembly","sequence","class","drug","score","level","mutations"]
    
    for opt, arg in opts:
        if opt in ('-h', '--help'):
            usage()
            sys.exit()
        elif opt in ('-i', '--input'):
            input = arg
        elif opt in ('-o', '--output'):
            output = arg
        else:
            usage()
            sys.exit()

    if not os.path.exists(input):
        print "Input file doesn't exist"
        sys.exit()
        
    if not os.path.exists(os.path.dirname(output)):
        try :
            os.makedirs(os.path.dirname(output))
        except OSError as exc:
            if exc.errno != errno.EEXIST:
                print "Missing or corrupted output dir"
                raise
    completedSRRID = {}
    cnt_records = 0
    if os.path.exists(output):
        with open(output) as output_fh:
            input_reader = csv.DictReader(output_fh)
            for row in input_reader:
                cnt_records += 1
                if row["srrID"] not in completedSRRID:
                    completedSRRID[row["srrID"]] = {}
                completedSRRID[row["srrID"]][row["assembly"]] = True
            
    cnt_inputs = 0    
    with open(input) as input_fh:
        input_reader = csv.DictReader(input_fh)
        for row in input_reader:
            cnt_inputs += 1
    
    current_input = 0
    with open(input) as input_fh, open(output, 'a') as output_fh:
        input_reader = csv.DictReader(input_fh)
        output_writer = csv.DictWriter(output_fh, header)
        
        if os.stat(output).st_size == 0:
            output_writer.writeheader()
        
        for row in input_reader:
            
            if not row["status log"].lower() == "completed":
                continue
            
            if row["srrID"] in completedSRRID and "contig" in completedSRRID[row["srrID"]]:
                print "Already processed contig of " + row["srrID"]
            elif "Contig path" in row:
                res = runHIVDBanalysis(row["Contig path"])
                cnt_records = printHIVDBresults(res, row["srrID"], "contig", cnt_records, output_writer)
                
            if row["srrID"] in completedSRRID and "global" in completedSRRID[row["srrID"]]:
                print "Already processed global of " + row["srrID"]
            elif "Global path" in row:
                res = runHIVDBanalysis(row["Global path"])
                cnt_records = printHIVDBresults(res, row["srrID"], "global", cnt_records, output_writer)
            
            current_input += 1
            print "Completed " + str(current_input) + " out of " + str(cnt_inputs)
            time.sleep(random.uniform(0.01,0.5))


if __name__ == '__main__':
    main()
