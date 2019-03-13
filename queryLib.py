#!/usr/bin/env python

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

def walkOneDir(workDir, fileExt):
    """Generates a list of all the file under workDir with fileExt as extention"""
    fileList = []
    filePattern = "*" + '.' + fileExt
    for pathName, dirName, fileNameList in os.walk(workDir):
        for fileName in fileNameList:
            if fnmatch.fnmatch(fileName, filePattern):
                fileList.append(os.path.abspath(os.path.join(pathName, fileName)))
    return fileList

def paraWalkDirBad(DirList, fileExt):
    """Walk directories to generate list of files parallelly"""
    procToRun = []
    pipeList = []
    fileList = []
    for subDir in DirList:
        recvEnd, sendEnd = multiprocessing.Pipe(False)
        procToRun.append(multiprocessing.Process(target=walkOneDir, args=(subDir, fileExt, sendEnd)))
        pipeList.append(recvEnd)
    for proc in procToRun:
        proc.start()
    for proc in procToRun:
        proc.join()
    for returnList in pipeList:
        fileList += returnList.recv()
    fileList = sorted(set(fileList))
    return fileList

def paraWalkDir(dirList, fileExt, threadNo):
    """Walk directories to generate list of files parallelly"""
    pool = multiprocessing.Pool(processes = threadNo)
    pipeList = []
    fileList = []
    for subDir in dirList:
        pipeList.append(pool.apply_async(walkOneDir, (subDir, fileExt)))
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
    for dirRec in recList:
        fwl.write("{0}\n".format(dirRec))
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
    recList = []
    sampList = []
    needUpdate = False
    for fileName in fileList:
        libName = os.path.basename(fileName)
        sampleName = libName.split('_')[0]
        sampList.append(sampleName)
        if not sampleName in fileDic.keys():
            fileDic[sampleName] = [fileName]
            needUpdate = True
        else:
            if not fileName in fileDic[sampleName]:
                fileDic[sampleName].append(fileName)
                needUpdate = True
        if needUpdate:
            recList.append(sampleName)
    fwd = open("{0}.del.log".format(dicFile), 'a')
    for sampleName in sampList:
        delLib = []
        oriLen = len(fileDic[sampleName])
        delLib = [libPath for libPath in fileDic[sampleName] if not os.access(libPath, os.F_OK)]
        fileDic[sampleName] = [libPath for libPath in fileDic[sampleName] if os.access(libPath, os.F_OK)]
        modLen = len(fileDic[sampleName])
        if modLen == 0:
            tmpKey = fileDic.pop(sampleName)
            needUpdate = True
        if oriLen != modLen:
            needUpdate = True
        if len(delLib) >= 1:
            fwd.write("{0}\n{1}\n".format(sampleName, "\n".join(delLib)))
    fwd.close()
    if needUpdate:
        fwl = open("{0}.rec.log".format(dicFile), 'a')
        for dirRec in recList:
            fwl.write("{0}\n".format(dirRec))
        fwl.close()
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

def regularyUpdate(libPath, dirLog, timeInterval, fileExt, threadNo, dicFile):
    """Regulary update database every timeInterval seconds"""
    updateLog = "{0}.upd.log".format(dicFile)
    while True:
        time.sleep(timeInterval)
        pathList = []
        nowTime = time.time()
        for dirName in os.listdir(libPath):
            pdir = os.path.join(libPath, dirName)
            if os.sccess(pdir, os.R_OK) and abs(nowTime, os.path.getmtime(pdir)) < 2 * timeInterval:
                for subDirName in os.listdir(pdir):
                    pSubDir = os.path.join(pdir, subDirName)
                    if os.access(pdir, os.R_OK) and abs(nowTime, os.path.getmtime(pSubDir)) < 2 * timeInterval:
                        pathList.append(pSubDir)
        pathList = sorted(set(pathList))
        if len(pathList) > 0:
            fw = open(updateLog, 'a')
            for pathName in pathList:
                fw.write("update database @ {0} from {1}\n".format(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), pathName))
            fw.close()
            fileList = paraWalkDir(pathList, fileExt, threadNo)
            updateDic(fileList, dicFile)
    return None

def forceUpdate(libPath, fileExt, threadNo, dicFile):
    """Force update specific lirary path"""
    dir3List = iniFileDir(libPath)
    fileList = paraWalkDir(dir3List, fileExt, threadNo)
    updateDic(fileList, dicFile)
    return None

def iniFileDir(fqPath):
    """ get directory path up to level 3 depth"""
    dirSet = set()
    dir1list = []
    for dirName in os.listdir(fqPath):
        pdir = os.path.join(fqPath, dirName)
        if os.path.isdir(pdir) and (not os.path.islink(pdir)) and os.access(pdir, os.R_OK):
            dir1list.append(pdir)
   
    dir2list = []
    for pdir in dir1list:
        for dir2 in os.listdir(pdir):
            pdir2 = os.path.join(pdir, dir2)
            if os.path.isdir(pdir2) and (not os.path.islink(pdir2)) and os.access(pdir2, os.R_OK):
                dir2list.append(pdir2)
            if (not os.path.isdir(pdir2)) and pdir2.endswith("clean.fastq.gz"):
                dirSet.add(pdir)
    
    dir3list = []
    for pdir in dir2list:
        for dir3 in os.listdir(pdir):
            pdir3 = os.path.join(pdir, dir3)
            if os.path.isdir(pdir3) and (not os.path.islink(pdir3)) and os.access(pdir3, os.R_OK):
                dirSet.add(pdir3)
            if (not os.path.isdir(pdir3)) and pdir3.endswith("clean.fastq.gz"):
                dirSet.add(pdir)
    return list(dirSet)

def main():
    """main function to accept arguments"""
    fqPath = "/share/seq_dir/ngs/"
    dicFile = "/home/wulj/database/queryLibDB/queryLibDB.json"
    fileExt = "clean.fastq.gz"
    threadNo = 10
    timeInterval = 3600
    userName = getpass.getuser()
    
    parser = argparse.ArgumentParser(usage=__doc__, formatter_class=SmartFormatter)
    parser.add_argument('-i', '--splist', help = textwrap.dedent("sample list, SampleName~FlowCell~..."), required = False)
    parser.add_argument('-o', '--outdir', help = textwrap.dedent("R|Output directory to put query log \ndefault: {0}".format(os.getcwd())), required = False)
    parser.add_argument('-t', '--jobtype', help = textwrap.dedent("R|Job type(q for query, c for create, u for update, s for summary) \ndefault: q"), required = False)
    parser.add_argument('-d', '--database', help = textwrap.dedent("R|Database file to create, update or query \ndefault: {0}".format(dicFile)), required = False)
    parser.add_argument('-f', '--fqpath', help = textwrap.dedent("R|Fastq path to create/update database \ndefault: {0}".format(fqPath)), required = False)
    parser.add_argument('-c', '--force', help = textwrap.dedent("R|Force update database"), required = False)
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
        dirList = iniFileDir(fqPath)
        fw = open("{0}.dir.log".format(dicFile), 'w')
        for dirName in dirList:
            fw.write("{0}\n".format(dirName))
        fw.close()
        fileList = paraWalkDir(dirList, fileExt, threadNo)
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
            forceUpdate(fqPath, dirLog, fileExt, threadNo, dicFile)
        else:
            regularyUpdate(fqPath, dirLog, timeInterval, fileExt, threadNo, dicFile)
    else:
        if args.database != None:
            dicFile = os.path.abspath(args.database)
            if not os.path.isfile(dicFile):
                print("Error, database file {0} is invalid!".format(dicFile))
                sys.exit(1)
        summaryDB(dicFile)

if __name__ == '__main__':
    main() 
