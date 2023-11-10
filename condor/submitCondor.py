import numpy as np 
import os
import glob
import argparse

d_args_debug = {
 200000:"mc.MGPy8EG_SVJSChan2j_750_2.py"
,200001:"mc.MGPy8EG_SVJSChan2j_750_8.py"
,200002:"mc.MGPy8EG_SVJSChan2j_2000_2.py"
,200003:"mc.MGPy8EG_SVJSChan2j_2000_8.py"
,200004:"mc.MGPy8EG_SVJSChan2j_5000_2.py"
,200005:"mc.MGPy8EG_SVJSChan2j_5000_8.py"
,200013:"mc.MGPy8EG_SVJSChan2j_1000_2.py"
,200014:"mc.MGPy8EG_SVJSChan2j_1000_2.py"
,200015:"mc.MGPy8EG_SVJSChan2j_1000_2.py"
,200018:"mc.MGPy8EG_SVJSChan2j_750_2.py"
,200019:"mc.MGPy8EG_SVJSChan2j_750_8.py"
,200020:"mc.MGPy8EG_SVJSChan2j_4000_2.py"
,200021:"mc.MGPy8EG_SVJSChan2j_4000_8.py"
,200022:"mc.MGPy8EG_SVJSChan2j_750_2.py"
,200023:"mc.MGPy8EG_SVJSChan2j_750_8.py"
,200024:"mc.MGPy8EG_SVJSChan2j_4000_2.py"
,200025:"mc.MGPy8EG_SVJSChan2j_4000_8.py"
,200026:"mc.MGPy8EG_SVJSChan2j_750_2.py"
,200027:"mc.MGPy8EG_SVJSChan2j_750_8.py"
,200028:"mc.MGPy8EG_SVJSChan2j_4000_2.py"
,200029:"mc.MGPy8EG_SVJSChan2j_4000_8.py"
}

d_args ={
 100000:"mc.MGPy8EG_SVJSChan2j_500_2.py"
,100001:"mc.MGPy8EG_SVJSChan2j_500_4.py"
,100002:"mc.MGPy8EG_SVJSChan2j_500_6.py"
,100003:"mc.MGPy8EG_SVJSChan2j_500_8.py"
,100004:"mc.MGPy8EG_SVJSChan2j_750_2.py"
,100005:"mc.MGPy8EG_SVJSChan2j_750_4.py"
,100006:"mc.MGPy8EG_SVJSChan2j_750_6.py"
,100007:"mc.MGPy8EG_SVJSChan2j_750_8.py"
,100008:"mc.MGPy8EG_SVJSChan2j_1000_2.py"
,100009:"mc.MGPy8EG_SVJSChan2j_1000_4.py"
,100010:"mc.MGPy8EG_SVJSChan2j_1000_6.py"
,100011:"mc.MGPy8EG_SVJSChan2j_1000_8.py"
,100012:"mc.MGPy8EG_SVJSChan2j_1250_2.py"
,100013:"mc.MGPy8EG_SVJSChan2j_1250_4.py"
,100014:"mc.MGPy8EG_SVJSChan2j_1250_6.py"
,100015:"mc.MGPy8EG_SVJSChan2j_1250_8.py"
,100016:"mc.MGPy8EG_SVJSChan2j_1500_2.py"
,100017:"mc.MGPy8EG_SVJSChan2j_1500_4.py"
,100018:"mc.MGPy8EG_SVJSChan2j_1500_6.py"
,100019:"mc.MGPy8EG_SVJSChan2j_1500_8.py"
,100020:"mc.MGPy8EG_SVJSChan2j_2000_2.py"
,100021:"mc.MGPy8EG_SVJSChan2j_2000_4.py"
,100022:"mc.MGPy8EG_SVJSChan2j_2000_6.py"
,100023:"mc.MGPy8EG_SVJSChan2j_2000_8.py"
,100024:"mc.MGPy8EG_SVJSChan2j_2500_2.py"
,100025:"mc.MGPy8EG_SVJSChan2j_2500_4.py"
,100026:"mc.MGPy8EG_SVJSChan2j_2500_6.py"
,100027:"mc.MGPy8EG_SVJSChan2j_2500_8.py"
,100028:"mc.MGPy8EG_SVJSChan2j_3000_2.py"
,100029:"mc.MGPy8EG_SVJSChan2j_3000_4.py"
,100030:"mc.MGPy8EG_SVJSChan2j_3000_6.py"
,100031:"mc.MGPy8EG_SVJSChan2j_3000_8.py"
,100032:"mc.MGPy8EG_SVJSChan2j_3500_2.py"
,100033:"mc.MGPy8EG_SVJSChan2j_3500_4.py"
,100034:"mc.MGPy8EG_SVJSChan2j_3500_6.py"
,100035:"mc.MGPy8EG_SVJSChan2j_3500_8.py"
,100036:"mc.MGPy8EG_SVJSChan2j_4000_2.py"
,100037:"mc.MGPy8EG_SVJSChan2j_4000_4.py"
,100038:"mc.MGPy8EG_SVJSChan2j_4000_6.py"
,100039:"mc.MGPy8EG_SVJSChan2j_4000_8.py"
,100040:"mc.MGPy8EG_SVJSChan2j_5000_2.py"
,100041:"mc.MGPy8EG_SVJSChan2j_5000_4.py"
,100042:"mc.MGPy8EG_SVJSChan2j_5000_6.py"
,100043:"mc.MGPy8EG_SVJSChan2j_5000_8.py"
,100044:"mc.MGPy8EG_SVJSChan2j_6000_2.py"
,100045:"mc.MGPy8EG_SVJSChan2j_6000_4.py"
,100046:"mc.MGPy8EG_SVJSChan2j_6000_6.py"
,100047:"mc.MGPy8EG_SVJSChan2j_6000_8.py"
,100048:"mc.MGPy8EG_SVJSChan2j_500_0.py"
,100049:"mc.MGPy8EG_SVJSChan2j_500_10.py"
,100050:"mc.MGPy8EG_SVJSChan2j_1500_0.py"
,100051:"mc.MGPy8EG_SVJSChan2j_1500_10.py"
,100052:"mc.MGPy8EG_SVJSChan2j_2000_0.py"
,100053:"mc.MGPy8EG_SVJSChan2j_2000_10.py"
,100054:"mc.MGPy8EG_SVJSChan2j_3000_0.py"
,100055:"mc.MGPy8EG_SVJSChan2j_3000_10.py"
,100056:"mc.MGPy8EG_SVJSChan2j_5000_0.py"
,100057:"mc.MGPy8EG_SVJSChan2j_5000_10.py"
,515519:"mc.MGPy8EG_SVJSChan2j_5000_2.py"
,515520:"mc.MGPy8EG_SVJSChan2j_5000_4.py"
,515521:"mc.MGPy8EG_SVJSChan2j_5000_6.py"
,515522:"mc.MGPy8EG_SVJSChan2j_5000_8.py"
,515523:"mc.MGPy8EG_SVJSChan2j_6000_2.py"
,515524:"mc.MGPy8EG_SVJSChan2j_6000_4.py"
,515525:"mc.MGPy8EG_SVJSChan2j_6000_6.py"
,515526:"mc.MGPy8EG_SVJSChan2j_6000_8.py"
,515593:"mc.MGPy8EG_SVJSChan2j_2500_2.py"
,515596:"mc.MGPy8EG_SVJSChan2j_2500_8.py"
,515603:"mc.MGPy8EG_SVJSChan2j_2500_2.py"
,515606:"mc.MGPy8EG_SVJSChan2j_2500_8.py"
,515613:"mc.MGPy8EG_SVJSChan2j_2500_2.py"
,515616:"mc.MGPy8EG_SVJSChan2j_2500_8.py"
}

#-------------------------------------------------------------------------
if __name__ == "__main__":
   parser = argparse.ArgumentParser()

   n_events = 10000 #10000

   #for i in [603,606,613,616]:  
   for i in range(0,58):  #88
     numName = '1000'
     if i < 10: numName += '0'+str(i)
     else: numName += str(i)
     #numName ='999999'

     # SUBMIT HERE
     print("SUBMITTING: input ", numName)
     args = open("args.txt","write")
     os.system("echo '"+numName+" "+str(n_events)+"' >>  args.txt")
     args.close()
     open("submit.sub","write")

     os.system("echo '#!/bin/bash' >> submit.sub")
     os.system("echo 'executable            = condor_evnt.sh' >> submit.sub") ##expand here to PC 
     os.system("echo 'output                = logs.$(ClusterId).$(ProcId).out' >> submit.sub")
     os.system("echo 'error                 = logs.$(ClusterId).$(ProcId).err' >> submit.sub")
     os.system("echo 'log                   = logs.$(ClusterId).log' >> submit.sub")
     os.system("echo 'universe         = vanilla' >> submit.sub")
     os.system("echo 'getenv           = True' >> submit.sub")
     os.system("echo 'Rank            = Mips' >> submit.sub")
     os.system("echo '+JobFlavour	= \"workday\"' >> submit.sub")
     os.system("echo '' >> submit.sub")
     os.system("echo 'should_transfer_files = YES' >> submit.sub")
     os.system("echo 'when_to_transfer_output = ON_EXIT' >> submit.sub")
     #os.system("echo 'initialdir = /afs/cern.ch/work/e/ebusch/public/SVJ/signal_request/condor/"+numName+"' >> submit.sub")
     os.system("echo 'initialdir = /afs/cern.ch/work/e/ebusch/public/SVJ/signal_request/condor/"+numName+"' >> submit.sub")
     #os.system("echo 'sampledir = /nevis/xenia/data/users/jgonski/xbb/Xbb_merged_samples/0121_PCJKDL1r' >> submit.sub")
     os.system("echo 'workdir = /afs/cern.ch/work/e/ebusch/public/SVJ/signal_request' >> submit.sub")
     os.system("echo 'transfer_input_files = $(workdir)/condor/condor_evnt.sh, $(workdir)/"+numName+"/"+d_args[int(numName)]+"' >> submit.sub")
     os.system("echo 'transfer_output_files = log.generate, DAOD_TRUTH1."+numName+".root' >> submit.sub")
     os.system("echo 'queue arguments from args.txt' >> submit.sub")

     os.system("condor_submit submit.sub")
     #time.sleep(.2)
     
     #open('submit.sub', 'w').close()

   print("DONE SUBMITTING... ")
