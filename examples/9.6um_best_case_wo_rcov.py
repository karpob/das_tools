#!/usr/bin/env python3
import os,sys,shutil,glob
# force matplotlib to use Agg backend to avoid broken X display stuff.
import matplotlib
matplotlib.use('Agg')
# find the path to this script.
thisDir = os.path.dirname(os.path.abspath(__file__))
# find the path the the python_tools path.
python_tools_path = os.path.join( os.path.dirname(thisDir),'python_tools')
#add it to the path, so we can import the library.
sys.path.insert(0,python_tools_path)
#import the 1d and 2d stat plotters
from stat_plotter import plotter_2d as plt2d
from stat_plotter import plotter_1d as plt1d

if __name__ == "__main__":
    startdate = 2018080100
    enddate = 2018093000
    exp_ids = ['bmk_o4_nf.21z']
    control = 'bmk_ctl_joff.21z'
    exp_names = ['9.6 $\mu$m']
    control_name = 'Control'
    domains = ['global','n.hem','s.hem','tropics']
    graphicOutput = '.png'
    confidence = 0.95
    f2move = '*'+graphicOutput
    storDir = os.path.join(thisDir, '9.6um_best_case_wo_rcov') 

    a = plt1d()
    a. plotter_1d(  startdate=startdate,\
                    enddate=enddate,\
                    exp_ids=exp_ids,\
                    control=control,\
                    exp_names=exp_names,\
                    control_name='Control',\
                    variables=['h','t','u','v','q'],\
                    domains=domains,\
                    confidence=confidence,\
                    verif1='gmao',\
                    verif2='gmao',\
                    maxlev=1000.,\
                    minlev=100.,\
                    levs=None,\
                    statistic='cor',\
                    graphicOutput=graphicOutput)

    outputDir = os.path.join( storDir, '1d_cor' )  
    if( not os.path.isdir(outputDir) ): os.makedirs( outputDir )
    files = glob.glob(os.path.join(thisDir,f2move))
    for f in files:
        shutil.move(f, os.path.join(outputDir,os.path.basename(f)) )
    a. plotter_1d(  startdate=startdate,\
                    enddate=enddate,\
                    exp_ids=exp_ids,\
                    control=control,\
                    exp_names=exp_names,\
                    control_name='Control',\
                    variables=['h','t','u','v','q'],\
                    domains=domains,\
                    confidence=confidence,\
                    verif1='gmao',\
                    verif2='gmao',\
                    maxlev=1000.,\
                    minlev=100.,\
                    levs=None,\
                    statistic='rms',\
                    graphicOutput=graphicOutput)

    outputDir = os.path.join( storDir, '1d_rms' )  
    if( not os.path.isdir(outputDir) ): os.makedirs( outputDir )
    files = glob.glob(os.path.join(thisDir,f2move))
    for f in files:
        shutil.move(f, os.path.join(outputDir,os.path.basename(f)) )


    b = plt2d()
    for i,e in enumerate(exp_ids):
        b.plotter_2d(   startdate=startdate,\
                        enddate=enddate,\
                        expid1=e,\
                        expid2=control,\
                        expid1_name=exp_names[i],\
                        expid2_name='Control',\
                        variables=['h','t','u','v','q'],\
                        domains=domains,\
                        confidence=confidence,\
                        sigplot=True,\
                        verif1='gmao',\
                        verif2='gmao',\
                        dblevs=None,\
                        maxlev=1000.,\
                        minlev=100.,\
                        levs=None,\
                        statistic='cor',\
                        graphicOutput=graphicOutput)

        outputDir = os.path.join( storDir, '2d_cor' )
  
        if( not os.path.isdir(outputDir) ): os.makedirs( outputDir )
        files = glob.glob(os.path.join(thisDir,f2move))
        for f in files:
            shutil.move(f, os.path.join(outputDir,os.path.basename(f)) )

        b.plotter_2d(  startdate=startdate,\
                        enddate=enddate,\
                        expid1=e,\
                        expid2=control,\
                        expid1_name=exp_names[i],\
                        expid2_name='Control',\
                        variables=['h','t','u','v','q'],\
                        domains=domains,\
                        confidence=confidence,\
                        sigplot=True,\
                        verif1='gmao',\
                        verif2='gmao',\
                        dblevs=None,\
                        maxlev=1000.,\
                        minlev=100.,\
                        levs=None,\
                        statistic='rms',\
                        graphicOutput=graphicOutput)

        outputDir = os.path.join( storDir, '2d_rms' )  
        if( not os.path.isdir(outputDir) ): os.makedirs( outputDir )
        files = glob.glob(os.path.join(thisDir,f2move))
        for f in files:
            shutil.move(f, os.path.join( outputDir, os.path.basename(f) ) )
