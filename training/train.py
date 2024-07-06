from __future__ import division
from __future__ import print_function

import time
import argparse
import numpy as np
import os

prime=[1009,1013,1019,1021,1031,1033,1039]


parser = argparse.ArgumentParser()
parser.add_argument('--seed', type=int, default=42, help='Random seed.')
parser.add_argument('--iter', type=int, default=200,
                    help='Number of epochs to train.')


args = parser.parse_args()


np.random.seed(args.seed)

#adj, features, labels, idx_train, idx_val, idx_test = load_data()
#Dataset=Synthetic_data("./syn_data/")

name="random-32-32-10"
map="instances/mapf-map/%s.map"%(name)
scen="instances/scen-random/%s-random"%(name)
num_of_agent=200
instance_set=[(i,prime[j]) for i in range(1,22) if i!=9 for j in range(0,1)]
train_set=[(i,prime[j])  for i in range(1,18) if i!=9 for j in range(0,1)]
val_set=[(i,prime[j]) for i in range(18,22) for j in range(0,1)]

time_limit=60
neighborsize=8
num_of_sample=20
solutionFile="feat/iter%d/solution_%d_%d.txt"
featFile="feat/iter%d/feature_%d_%d.txt"

os.system("mkdir feat")
os.system("mkdir feat/iter0")
os.system("mkdir exp")


def runLNS(iter):
    #return
    #if iter<=0: return
    if iter==0:
        initSol="NONE"
    else:
        initSol=""
    __=0
    for i,rand in instance_set:
        if initSol!="NONE": 
            initSol=solutionFile%(iter-1,i,rand)
        outputSol=solutionFile%(iter,i,rand)
        featfile=featFile%(iter,i,rand)
        #fixed neighborhood size: runLNS="./lns -m %s -a %s-%d.scen -k %d -t %d --neighborSize=%d --maxIterations=2 --initAlgo=PPS --destoryStrategy=Adaptive --replanAlgo=PP -o exp/%s-%d.csv --stats=exp/%s-%d_stats.csv "%(map,scen,i,num_of_agent,time_limit,neighborsize,name,iter,name,iter)
        #runCBSH2="./CBSH2/CBSH2 -m %s -a %s-%d.scen -t 60 -s 1 -h WDG -k %d --rand %d --featureFile %s"%(map,scen,i,num_of_agent,rand,featfile+"GraphDG")
        runLNS="./lns -m %s -a %s-%d.scen -k %d -t %d --lbns=5 --ubns=16 --maxIterations=2 --initAlgo=PPS --destoryStrategy=Adaptive --replanAlgo=PP -o exp/%s-%d.csv --stats=exp/%s-%d_stats.csv "%(map,scen,i,num_of_agent,time_limit,name,i,name,i)
        tp="--featureFile=%s --samples=%d --rand=%d --initSol=%s --outputSol=%s"%(featfile,num_of_sample,rand,initSol,outputSol)
        runLNS+=tp
        #print(runCBSH2)
        if __==0: print(runLNS)
        #os.system(runCBSH2)
        os.system(runLNS)
        __+=1
        #if __==5: return
        if __%25==0: print("Finish collecting %d instances"%(__))
    print("Finish data collection")



global QID
QID=1
#global outputfile
#outputfile="training/process_features/val_data.txt"
#instance_set=val_set

def process_feature_data(path,instance_set,output_file,validate,model="NA"):
    '''
    with open(path+"Graph","r") as f:
        line=[x for x in f]
        source=eval(line[0])
        target=eval(line[1])
        weight=eval(line[2])
    '''
    Rank=dict()
    Neighbor=dict()
    cnt_improved=0
    cnt_can_improve=0
    improved_value=0
    can_improve_value=0
    for a,b in instance_set:
        data_list=[]
        val=[]
        norm=[0]*1000
        max_score=0
        score_list=[]
        for id in range(num_of_sample):
            featfile=path+"feature_%d_%d.txt%d"%(a,b,id)
            if os.stat(featfile).st_size==0:  print("no feature!!!")
            x=[]
            y=[]
            neighbor=[]
            data=[]
            feat=[[[] for __ in range(100)] for _ in range(2)]
            feat_len=0
            with open(featfile,"r") as f:
                line=[x for x in f]
                agent=eval(line[0])
                line=line[1:]
                for i in range(agent):
                    node_feat=eval(line[i])
                    x.append(node_feat)
                    neighbor.append(node_feat[0])
                    k=int(node_feat[0])
                    feat_len=len(node_feat)-1
                    for j in range(1,len(node_feat)):
                        feat[k][j-1].append(node_feat[j])
                

                y=eval(line[agent])
                
                #score=y[0]*1.0
                score=y[1]*1.0
                data.append(sum(neighbor))
                for i in [0,1]:
                    for j in range(feat_len):
                        data.append(min(feat[i][j]))
                        data.append(max(feat[i][j]))
                        data.append(sum(feat[i][j]))
                        data.append(sum(feat[i][j])/len(feat[i][j]))
                for i in range(len(data)):
                    norm[i]=max(norm[i],abs(data[i]))
                data_list.append(data)
                val.append([score,0,neighbor])
                score_list.append(score)
                max_score=max(max_score,score)

        score_list.sort(reverse=True)
        score_margin=max(max_score*0.6,score_list[int(num_of_sample*0.2)-1])
        score_margin2=score_list[int(num_of_sample*0.5)-1]
        global QID
        if validate: 
            with open(output_file,"w") as f:
                __adf=1
        with open(output_file,"a") as f:
            for id in range(num_of_sample):
                if not validate:
                    if val[id][0]>=score_margin:
                        f.write("2 qid:%d"%(QID))
                    elif val[id][0]>=score_margin2:
                        f.write("1 qid:%d"%(QID))
                    else:
                        f.write("0 qid:%d"%(QID))
                else:
                    f.write("%.4lf qid:%d"%(val[id][0],1))
                for i in range(len(data_list[id])):
                    if norm[i]<1e-5: norm[i]=1
                    f.write(" %d:%.4lf"%(i+1,data_list[id][i]/norm[i]))
                f.write("\n")
            QID+=1
        if not validate: 
            continue
        #print(a,b)
        #print("./training/svm_rank_classify %s %s training/prediction"%(output_file,model))
        os.system("./training/svm_rank_classify %s %s training/prediction"%(output_file,model))
        with open("training/prediction","r") as f:
            line=[eval(x) for x in f]
            max_pred=0
            for i in range(len(line)):
                if line[max_pred]<line[i]:
                    max_pred=i
            pred_rank=1
            if val[max_pred][0]>0: cnt_improved+=1
            if score_list[0]>0: cnt_can_improve+=1
            improved_value+=val[max_pred][0]
            can_improve_value+=score_list[0]
            #print(score_list,max_pred)
            for i in range(len(score_list)):
                if score_list[i]>val[max_pred][0]: pred_rank+=1
            Rank[(a,b)]=pred_rank
            Neighbor[(a,b)]=val[max_pred][2]
        
    return Rank,Neighbor,cnt_improved,cnt_can_improve,improved_value,can_improve_value

t_total = time.time()
Titer = args.iter


def validate():
    #with open("training/val_log.txt","w") as f:
    #    __WEAf=0
    for i in range(1,Titer+1):
        with open("training/val_log.txt","a") as f:
            f.write(str(i)+":\n")
        tot_rank=0
        cnt=0
        c1=0
        c2=0
        c3=0
        c4=0
        for j in range(1,(Titer+1)):
            rank,neighbor,cnt_imp,cnt_can,val_imp,val_can=process_feature_data("feat/iter"+str(j-1)+"/",val_set,"training/data/validate_temp.txt",True,"training/saved_model/model"+str(i))
            #rank,neighbor=process_feature_data("feat/iter"+str(j-1)+"/",val_set,"training/data/validate_temp.txt",True,"training/model1")
            c1+=cnt_imp
            c2+=cnt_can
            c3+=val_imp
            c4+=val_can
            for a,b in val_set:
                tot_rank+=rank[(a,b)]
                cnt+=1
            if j%20==0:
                with open("training/val_log.txt","a") as f:
                    f.write(str(j)+" "+str(tot_rank/cnt)+" improved cnt "+str(c1/c2)+" improved value "+str(c3/c4)+"\n")
        print(str(i)+" "+str(cnt)+" ave_rank:"+str(tot_rank/cnt)+"\n")
        #return
        
        with open("training/val_log.txt","a") as f:
            f.write(str(cnt)+" ave_rank:"+str(tot_rank/cnt)+"\n")


os.system("mkdir training/saved_model")
os.system("mkdir training/data")

model_list=[]
runLNS(0)

with open("training/trainlog.txt","w") as f:
    xasdf_=0
with open("training/data/training_data.txt","w") as f:
    xasdf_=0

Dataset=[]
Dataset_val=[]

ave_rank_val=[0]*1000
ave_rank_train=[0]*1000


for i in range(1,Titer+1):
    #if i==1:
        #runLNS(1)
        #continue
    os.system("mkdir feat/iter"+str(i))
    os.system("mkdir training/syn_data/iter"+str(i))
    os.system("mkdir training/syn_data/iter_val"+str(i))
    process_feature_data("feat/iter"+str(i-1)+"/",train_set,"training/data/training_data.txt",False)
    TrainCommand="./training/svm_rank_learn -c 10 -t 0 -e 0.01 -d 2 -s 1 -r 1 -l 2 training/data/training_data.txt training/saved_model/model"+str(i)
    os.system(TrainCommand)
    
    os.system("sleep 10")

    rank,neighbor,_1,_2,_3,_4=process_feature_data("feat/iter"+str(i-1)+"/",instance_set,"training/data/validate_temp.txt",True,"training/saved_model/model"+str(i))
    for a,b in instance_set:
        outputSol=solutionFile%(i-1,a,b)
        #print(a,rank[a,b])
        if a>=17: 
            ave_rank_val[i]+=rank[a,b]/4.
        else:
            ave_rank_train[i]+=rank[a,b]/16.
        with open(outputSol,"a+") as f:
            #neighbor=[]
            f.write("%d"%(sum(neighbor[a,b])))
            #print(outputSol,sum(pred))
            for k in range(len(neighbor[a,b])):
                if neighbor[a,b][k]==1:f.write(" %d"%(k))
            f.write("\n")
    log=[(a,rank[a,b]) for a,b in instance_set]
    print(log)
    print('Iter: {:03d}, ave_rank_train: {:.5f} , ave_rank_val: {:.5f}\n'.
            format(i, ave_rank_train[i],ave_rank_val[i]))
    with open("training/trainlog.txt","a") as f:
        f.write(str(log))
        f.write('Iter: {:03d}, ave_rank_train: {:.5f} , ave_rank_val: {:.5f}\n'.
            format(i, ave_rank_train[i],ave_rank_val[i]))
        
            #label = data.y.detach().cpu().numpy()
    #print('Iter: {:03d}, Acc_train: {:.5f} , Precision_train: {:.5f}, Recall_train: {:.5f}'.format(i, acc_train/cnt,precision_train/cnt,recall_train/cnt))
    '''
    with open("trainlog.txt","a") as f:
        f.write('Iter: {:03d}, ave_rank_train: {:.5f} , ave_rank_val: {:.5f}\n'.
            format(i, ave_rank_train[i],ave_rank_val[i]))
    '''
    '''
    t=len(Dataset)-10
    if t<0:t=0
    if len(Dataset)>1:
        acc_train=0
        precision_train=0
        recall_train=0
        cnt=0
        eval_loader = DataLoader(ConcatDataset(Dataset[t:-1]), batch_size=1)
        with torch.no_grad():
            for data in eval_loader:
                data = data.to(device)
                pred = model(data).detach().cpu()#.numpy()
                pred = torch.argmax(pred,dim=1)
                if data.y.shape[0]>0:
                    cnt+=1
                    acc_train+=(pred==data.y).float().mean().item()
                    precision.update((pred,data.y))
                    pre=precision.compute()
                    precision_train+= pre.item()
                    recall.update((pred,data.y))
                    recall_train+=recall.compute().item()
                pred = pred.numpy()            


    print('Iter: {:03d}, Acc_train: {:.5f} , Precision_train: {:.5f}, Recall_train: {:.5f}'.
            format(i, acc_train/cnt,precision_train/cnt,recall_train/cnt))
    with open("trainlog.txt","a") as f:
        f.write('Iter: {:03d}, Acc_train: {:.5f} , Precision_train: {:.5f}, Recall_train: {:.5f}\n'.
            format(i, acc_train/cnt,precision_train/cnt,recall_train/cnt))
    '''

    #torch.save(model, "training/saved_model/iter%d.pt"%(i))
    #model_list.append("training/saved_model/iter%d"%(i))
    #runLNS(i)
    runLNS(i)
                

validate()
