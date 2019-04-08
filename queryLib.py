#!/usr/bin/env python3

import os
import re
import sys
import json
import time
import errno
import shutil
import fnmatch
import getpass
import argparse
import datetime
import textwrap
import multiprocessing

__doc__ = "script to build, update and query fastq path database"

class SmartFormatter(argparse.HelpFormatter):
    """Define SmartFormatter to display help info in multiple lines"""
    def _split_lines(self, text, width):
        if text.startswith('R|'):
            return text[2:].splitlines()  
        return argparse.HelpFormatter._split_lines(self, text, width)

def walkOneDir(workDir, fileExtPat):
    """Generates a list of all the files under workDir with fileExtPat as extention"""
    fileList = []
    for pathName, dirName, fileNameList in os.walk(workDir):
        for fileName in fileNameList:
            if re.match(fileExtPat, fileName):
                fileList.append(os.path.abspath(os.path.join(pathName, fileName)))
    return fileList

def paraWalkDir(dirList, fileExtPat, threadNo):
    """Walk directories to generate list of files parallelly"""
    pool = multiprocessing.Pool(processes = threadNo)
    pipeList = []
    fileList = []
    for subDir in dirList:
        pipeList.append(pool.apply_async(walkOneDir, (subDir, fileExtPat)))
    pool.close()
    pool.join()
    for returnList in pipeList:
        fileList += returnList.get()
    fileList = sorted(set(fileList))
    return fileList

def makeFileDic(fileList, dicFile):
    """Make dictinoary of {SampleName:[libraryPathA...libraryPathX]}"""
    fileDic = {}
    recList = []
    for fileName in fileList:
        libName = os.path.basename(fileName)
        sampleName = libName.split('_')[0]
        if not sampleName in fileDic.keys():
            fileDic[sampleName] = [fileName]
        else:
            fileDic[sampleName].append(fileName)
        recList.append(fileName)
    recList = sorted(set(recList))
    fwl = open("{0}.rec.log".format(dicFile), 'w')
    fwl.write("{0}\n".format("\n".join(recList)))
    fwl.close()
    fw = open(dicFile, 'w')
    json.dump(fileDic, fw)
    fw.close()
    return None

def updateDic(fileList, dicFile):
    """Update dictionary json file"""
    fr = open(dicFile, 'r')
    fileDic = json.load(fr)
    fr.close()
    dicKeys = fileDic.keys()
    # check for new values needed to inserted
    recList = []
    for fileName in fileList:
        needUpdate = False
        libName = os.path.basename(fileName)
        sampleName = libName.split('_')[0]
        if not sampleName in dicKeys:
            fileDic[sampleName] = [fileName]
            needUpdate = True
        else:
            if not fileName in fileDic[sampleName]:
                fileDic[sampleName].append(fileName)
                needUpdate = True
        if needUpdate:
            recList.append(fileName)
    if len(recList) >= 1:
        fwl = open("{0}.rec.log".format(dicFile), 'a')
        fwl.write("{0}\n".format("\n".join(recList)))
        fwl.close();
    # check for old values needed to be removed
    delLib = []
    dropkey = []
    for sampleName in dicKeys:
        oriLen = len(fileDic[sampleName])
        delLib.extend([libPath for libPath in fileDic[sampleName] if not os.access(libPath, os.F_OK)])
        fileDic[sampleName] = [libPath for libPath in fileDic[sampleName] if os.access(libPath, os.F_OK)]
        modLen = len(fileDic[sampleName])
        if modLen == 0:
            dropkey.append(sampleName);
    for k in dropkey:
        fileDic.pop(k)
    if len(delLib) >= 1:
        fwd = open("{0}.del.log".format(dicFile), 'a')
        fwd.write("{0}\n{1}\n".format(sampleName, "\n".join(delLib)))
        fwd.close()
    # update json database file
    if len(recList) > 0 or len(delLib) > 0:
        newDicFile = "{0}.new".format(dicFile)
        fw = open(newDicFile, 'w')
        json.dump(fileDic, fw)
        fw.close()
        shutil.move(newDicFile, dicFile)
    return None
    
def queryFileDic(dicFile, sfDic, outPath, patList):
    """Query library dictionary file to get libraries"""
    fr = open(dicFile, 'r')
    fileDic = json.load(fr)
    fr.close()

    failLog = open(os.path.join(outPath, 'queryFail.log'), 'w')
    succLog = open(os.path.join(outPath, 'querySucc.log'), 'w')

    dicKeys = fileDic.keys()
    for sampleName in sfDic.keys():
        if sampleName in dicKeys:
            libGet = False
            matchList = []
            for libPath in fileDic[sampleName]:
                pattMatch = True
                finaPatList = [e for e in patList]
                for flowCell in sfDic[sampleName]:
                    finaPatList.append(".*" + flowCell + ".*")
                for libPatt in finaPatList:
                    if not re.match(libPatt, libPath):
                        pattMatch = False
                        break
                if pattMatch:
                    libGet = True
                    matchList.append(libPath)
            if libGet:
                succLog.write("{0}\t{1}\n".format(sampleName, "\t".join(matchList)))
            else:
                failLog.write("{0} not found in {1} to match your pattern!\n".format(sampleName, dicFile))   
        else:
            failLog.write("{0} not found in {1}!\n".format(sampleName, dicFile))
    return None

def parseSpList(sampleList):
    """extract SampleName from sampleList"""
    fr = open(sampleList, 'r')
    fr.readline()
    sfdic = {} #dictinary of SampleName : [FlowCell]
    for line in fr:
        tmpList = line.rstrip().split('\t')
        SampleName = tmpList[0]
        FlowCell = tmpList[1]
        if SampleName in sfdic.keys():
            sfdic[SampleName].add(FlowCell)
        else:
            sfdic[SampleName] = set()
            sfdic[SampleName].add(FlowCell)
    return sfdic

def regularyUpdate(libPath, blackList, dirLog, timeInterval, fileExtPat, threadNo, dicFile):
    """Regulary update database every timeInterval seconds"""
    while True:
        time.sleep(3600)
        pathList = getModFileDir(libPath, blackList, timeInterval, fileExtPat)
        if len(pathList) > 0:
            updateLog = "{0}.mup.log".format(dicFile)
            fw = open(updateLog, 'w')
            fw.write("paths may needed to be updated @ {0}:\n".format(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
            fw.write("{0}\n".format("\n".join(pathList)))
            fw.close()
            fileList = paraWalkDir(pathList, fileExtPat, threadNo)
            updateDic(fileList, dicFile)
    return None

def forceUpdate(libPath, blackList, fileExtPat, threadNo, dicFile):
    """Force update specific lirary path"""
    dir3List = iniFileDir(libPath, blackList, fileExtPat)
    fileList = paraWalkDir(dir3List, fileExtPat, threadNo)
    updateDic(fileList, dicFile)
    return None

def iniFileDir(fqPath, blackList, fileExtPat):
    """ get directory path up to level 3 depth"""
    dirSet = set()
    dir1list = []
    if not os.access(fqPath, os.F_OK):
        return dirSet
    for dirName in os.listdir(fqPath):
        pdir = os.path.abspath(os.path.join(fqPath, dirName))
        if (not pdir in blackList) and os.path.isdir(pdir) and (not os.path.islink(pdir)) and os.access(pdir, os.R_OK):
            dir1list.append(pdir)
   
    dir2list = []
    for pdir in dir1list:
        if not os.access(pdir, os.F_OK):
            continue
        for dir2 in os.listdir(pdir):
            pdir2 = os.path.join(pdir, dir2)
            if os.path.isdir(pdir2) and (not os.path.islink(pdir2)) and os.access(pdir2, os.R_OK):
                dir2list.append(pdir2)
            if (not os.path.isdir(pdir2)) and re.match(fileExtPat, pdir2):
                dirSet.add(pdir)
    
    dir3list = []
    for pdir in dir2list:
        if not os.access(pdir, os.F_OK):
            continue
        for dir3 in os.listdir(pdir):
            pdir3 = os.path.join(pdir, dir3)
            if os.path.isdir(pdir3) and (not os.path.islink(pdir3)) and os.access(pdir3, os.R_OK):
                dirSet.add(pdir3)
            if (not os.path.isdir(pdir3)) and re.match(fileExtPat, pdir3):
                dirSet.add(pdir)
    return list(dirSet)

def getModFileDir(fqPath, blackList, timeInterval, fileExtPat):
    """ get directory path up to level 3 depth"""
    dirSet = set()
    dir1list = []
    if not os.access(fqPath, os.F_OK):
        return dirSet
    for dirName in os.listdir(fqPath):
        pdir = os.path.abspath(os.path.join(fqPath, dirName))
        if (not pdir in blackList) and os.path.isdir(pdir) and (not os.path.islink(pdir)) and os.access(pdir, os.R_OK):
            dir1list.append(pdir)
   
    dir2list = []
    for pdir in dir1list:
        if not os.access(pdir, os.F_OK):
            continue
        for dir2 in os.listdir(pdir):
            pdir2 = os.path.join(pdir, dir2)
            if os.path.isdir(pdir2) and (not os.path.islink(pdir2)) and os.access(pdir2, os.R_OK):
                dir2list.append(pdir2)
            if (not os.path.isdir(pdir2)) and re.match(fileExtPat, pdir2) and abs(time.time() - os.path.getmtime(pdir2)) < timeInterval:
                dirSet.add(pdir)
    
    dir3list = []
    for pdir in dir2list:
        if not os.access(pdir, os.F_OK):
            continue
        for dir3 in os.listdir(pdir):
            pdir3 = os.path.join(pdir, dir3)
            if os.path.isdir(pdir3) and (not os.path.islink(pdir3)) and os.access(pdir3, os.R_OK) and abs(time.time() - os.path.getmtime(pdir3)) < timeInterval:
                dirSet.add(pdir3)
            if (not os.path.isdir(pdir3)) and re.match(fileExtPat, pdir3) and abs(time.time() - os.path.getmtime(pdir3)) < timeInterval:
                dirSet.add(pdir)
    return list(dirSet)

def main():
    """main function to accept arguments"""
    fqPath = "/share/seq_dir/ngs/"
    dicFile = "/share/work1/wulj/database/queryLibDB/queryLibDB.json"
    blackList = ["/share/seq_dir/ngs/Runs", "/share/seq_dir/ngs/qc"]

    fileExtPat = re.compile(".*R[12]+(.clean)?.fastq.gz$")
    threadNo = 10
    timeInterval = 36000
    userName = getpass.getuser()
    
    parser = argparse.ArgumentParser(usage=__doc__, formatter_class=SmartFormatter)
    parser.add_argument('-i', '--splist', help = textwrap.dedent("sample list, SampleName~FlowCell~..."), required = False)
    parser.add_argument('-o', '--outdir', help = textwrap.dedent("R|Output directory to put query log \ndefault: {0}".format(os.getcwd())), required = False)
    parser.add_argument('-t', '--jobtype', help = textwrap.dedent("R|Job type(q for query, c for create, u for update) \ndefault: q"), required = False)
    parser.add_argument('-d', '--database', help = textwrap.dedent("R|Database file to create, update or query \ndefault: {0}".format(dicFile)), required = False)
    parser.add_argument('-f', '--fqpath', help = textwrap.dedent("R|Fastq path to create/update database \ndefault: {0}".format(fqPath)), required = False)
    parser.add_argument('-c', '--force', help = textwrap.dedent("R|Force update database"), action = 'store_true')
    parser.add_argument('-g', '--libpatterns', help = textwrap.dedent("R|Library patterns to match in querying, seperated by comma \ndefault: *"), required = False)
    parser.add_argument('-p', '--thread', help = "R|Thread number used to create/update database \ndefault: {0}".format(threadNo), required = False)
    args = parser.parse_args()

    jobType = "q"
    if args.jobtype != None:
        jobType = args.jobtype.lower()
    if not jobType in ["c", "u", "q"]:
        print("Error, only c(create), u(update) or q(query) jobtype is allowed!")
        sys.exit(1)

    if jobType == "q":
        spList = ""
        outDir = ""
        if args.splist == None:
            print("Error, sample list is not provided for query!")
            sys.exit(1)
        else:
            spList = os.path.abspath(args.splist)
        if not os.path.isfile(spList):
            print("Error, sample list {0} doesn't exist!".format(spList))
            sys.exit(1)
        outDir = os.getcwd()
        if args.outdir != None:
            outDir = os.path.abspath(args.outdir)
            if not os.access(outDir, os.F_OK):
                os.mkdir(outDir)
        if args.database != None:
            dicFile = os.path.abspath(args.database)
            if not os.path.isfile(dicFile):
                print("Error, database file {0} is invalid!".format(dicFile))
                sys.exit(1)
        patList = []
        libPatt = ".*"
        if args.libpatterns != None:
            libPatt = args.libpatterns
        for rawPat in libPatt.split(","):
            patList.append(re.compile(".*{0}.*".format(rawPat)))
        sfDic = parseSpList(spList)
        queryFileDic(dicFile, sfDic, outDir, patList)
    elif jobType == "c":
        if args.fqpath != None:
            fqPath = os.path.abspath(args.fqpath)
            if not os.path.isdir(fqPath):
                print("Error, invalid fastq path:{0}".format(fqPath))
                sys.exit(1)
            else:
                print("fqpath:{0}".format(fqPath))
        
        if args.database != None:
            dicFile = os.path.abspath(args.database)
        dicDir = os.path.dirname(dicFile)
        if not os.access(dicDir, os.F_OK):
            os.mkdir(dicDir)
        if not os.access(dicFile, os.F_OK):
            open(dicFile, 'w').close()   
        print("database path:{0}".format(dicFile))
        if not os.access(dicFile, os.W_OK):
            print("{0}, you're not allowed to create this databae:{1}!".format(userName, dicFile))
            sys.exit(1)
        if args.thread != None:
            threadNo = int(args.thread)
        dirList = iniFileDir(fqPath, blackList, fileExtPat)
        fw = open("{0}.dir.log".format(dicFile), 'w')
        fw.write("{0}\n".format("\n".join(dirList)))
        fw.close()
        fileList = paraWalkDir(dirList, fileExtPat, threadNo)
        makeFileDic(fileList, dicFile)
    elif jobType == "u":
        if args.fqpath != None:
            fqPath = os.path.abspath(args.fqpath)
            if not os.path.isdir(fqPath):
                print("Error, invalid fastq path:{0}".format(fqPath))
                sys.exit(1)
        if args.thread != None:
            threadNo = int(args.thread)
        if args.database != None:
            dicFile = os.path.abspath(args.database)
            if not os.path.isfile(dicFile):
                print("Error, database file {0} is invalid!".format(dicFile))
                sys.exit(1)
        if not os.access(dicFile, os.W_OK):
            print("{0}, you're not allowed to update this databae:{1}!".format(userName, dicFile))
            sys.exit(1)
        dirLog = "{0}.dir.log".format(dicFile)
        if args.force:
            forceUpdate(fqPath, blackList, fileExtPat, threadNo, dicFile)
        else:
            regularyUpdate(fqPath, blackList, dirLog, timeInterval, fileExtPat, threadNo, dicFile)
    else:
        if args.database != None:
            dicFile = os.path.abspath(args.database)
            if not os.path.isfile(dicFile):
                print("Error, database file {0} is invalid!".format(dicFile))
                sys.exit(1)
        summaryDB(dicFile)

if __name__ == '__main__':
    main() 
