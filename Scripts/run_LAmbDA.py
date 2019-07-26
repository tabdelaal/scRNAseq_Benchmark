# -*- coding: utf-8 -*-
"""
Created on Thu May 23 13:51:15 2019

@author: Lieke
"""

import os 
import numpy as np
import pandas as pd
import time as tm
import rpy2.robjects as robjects
import tensorflow as tf
import math
import scipy.io as sio
import optunity as opt
from tensorflow.contrib.tensor_forest.python import tensor_forest
from tensorflow.python.ops import resources


def run_LAmbDA(DataPath, LabelsPath, CV_RDataPath, OutputDir, GeneOrderPath = "", NumGenes = 0):
    '''
    run LAmbDA classifier
    Wrapper script to run LAmbDA on a benchmark dataset with 5-fold cross validation,
    outputs lists of true and predicted cell labels as csv files, as well as computation time.
  
    Parameters
    ----------
    DataPath : Data file path (.csv), cells-genes matrix with cell unique barcodes 
    as row names and gene names as column names.
    LabelsPath : Cell population annotations file path (.csv).
    CV_RDataPath : Cross validation RData file path (.RData), obtained from Cross_Validation.R function.
    OutputDir : Output directory defining the path of the exported file.
    GeneOrderPath : Gene order file path (.csv) obtained from feature selection, 
    defining the genes order for each cross validation fold, default is NULL.
    NumGenes : Number of genes used in case of feature selection (integer), default is 0.
    '''
        
    # read the Rdata file
    robjects.r['load'](CV_RDataPath)

    nfolds = np.array(robjects.r['n_folds'], dtype = 'int')
    tokeep = np.array(robjects.r['Cells_to_Keep'], dtype = 'bool')
    col = np.array(robjects.r['col_Index'], dtype = 'int')
    col = col - 1 
    test_ind = np.array(robjects.r['Test_Idx'])
    train_ind = np.array(robjects.r['Train_Idx'])

    # read the data
    data = pd.read_csv(DataPath,index_col=0,sep=',')
    labels = pd.read_csv(LabelsPath, header=0,index_col=None, sep=',', usecols = col)
    
    labels = labels.iloc[tokeep]
    data = data.iloc[tokeep]
    
    # read the feature file
    if (NumGenes > 0):
        features = pd.read_csv(GeneOrderPath,header=0,index_col=None, sep=',')
    
    # folder with results
    os.chdir(OutputDir)
                
    tr_time=[]
    ts_time=[]
    truelab = np.zeros([len(labels),1],dtype = int)
    predlab = np.zeros([len(labels),1],dtype = int)
        
    for i in range(np.squeeze(nfolds)):
        global X, Y, Gnp, Dnp, train, test, prt, cv
        test_ind_i = np.array(test_ind[i], dtype = 'int') - 1
        train_ind_i = np.array(train_ind[i], dtype = 'int') - 1
                
        X = np.array(data) 
        if (NumGenes > 0):
            X = np.log2(X/10+1)
            feat_to_use = features.iloc[0:NumGenes,i]
            X = X[:,feat_to_use]
        else:
            X = np.log2(np.transpose(select_feats(np.transpose(X),0.5,80))/10+1)
    
        uniq = np.unique(labels)
        Y = np.zeros([len(labels),len(uniq)],int)
        
        for j in range(len(uniq)):
            Y[np.where(labels == uniq[j])[0],j] = 1
    
        Y = np.array(Y)
        
        Gnp = np.zeros([len(uniq),len(uniq)],int)
        np.fill_diagonal(Gnp,1)
        Gnp = np.array(Gnp)
        
        Dnp = np.ones([len(uniq),1],int)
        Dnp = np.array(Dnp)
        
        train_samp = int(np.floor(0.75*len(train_ind_i)))
        test_samp = len(train_ind_i) - train_samp
        perm = np.random.permutation(len(train_ind_i))
        train = perm[0:train_samp]
        test = perm[train_samp:test_samp+1]
        
        while(np.sum(np.sum(Y[train,:],0)<5)>0):
            perm = np.random.permutation(X.shape[0])
            train = perm[0:train_samp+1]
            test = perm[train_samp+1:train_samp+test_samp+1]
        
        cv = i
        optunity_it = 0
        prt = False
        opt_params = None
                    
        start=tm.time()
        opt_params, _, _ = opt.minimize(run_LAmbDA2,solver_name='sobol', gamma=[0.8,1.2], delta=[0.05,0.95], tau=[10.0,11.0], prc_cut=[20,50], bs_prc=[0.2,0.6], num_trees=[10,200], max_nodes=[100,1000], num_evals=50)
        tr_time.append(tm.time()-start)
        
        print("Finished training!")
        
        prt = True
        train = train_ind_i
        test = test_ind_i
        
        start=tm.time()
        err = run_LAmbDA2(opt_params['gamma'], opt_params['delta'], opt_params['tau'], opt_params['prc_cut'], opt_params['bs_prc'], opt_params['num_trees'], opt_params['max_nodes'])
        ts_time.append(tm.time()-start)
        
        tf.reset_default_graph();
        
        predfile = 'preds_cv' + str(cv) + '.mat'
        truefile = 'truth_cv' + str(cv) + '.mat'
        pred = sio.loadmat(predfile)
        truth = sio.loadmat(truefile)
        
        pred = pred['preds']
        truth = truth['labels']
        
        pred_ind = np.argmax(pred,axis=1)
        truth_ind = np.argmax(truth,axis=1)
        
        predlab[test_ind_i,0] = pred_ind
        truelab[test_ind_i,0] = truth_ind
            
                
    truelab = pd.DataFrame(truelab)
    predlab = pd.DataFrame(predlab)
        
    tr_time = pd.DataFrame(tr_time)
    ts_time = pd.DataFrame(ts_time)
        
    if (NumGenes == 0):  
        truelab.to_csv("LAmbDA_True_Labels.csv", index = False)
        predlab.to_csv("LAmbDA_Pred_Labels.csv", index = False)
        tr_time.to_csv("LAmbDA_Training_Time.csv", index = False)
        ts_time.to_csv("LAmbDA_Testing_Time.csv", index = False)
    else:
        truelab.to_csv("LAmbDA_" + str(NumGenes) + "_True_Labels.csv", index = False)
        predlab.to_csv("LAmbDA_" + str(NumGenes) + "_Pred_Labels.csv", index = False)
        tr_time.to_csv("LAmbDA_" + str(NumGenes) + "_Training_Time.csv", index = False)
        ts_time.to_csv("LAmbDA_" + str(NumGenes) + "_Testing_Time.csv", index = False)


##### Functions copied from LAmbDA's Github
def wt_cutoff(colnum,cutoff,Gtmp,gamma):
	rowsums = np.sum(Gtmp,axis=1);
	return(math.ceil(cutoff*(math.log((max(rowsums)/rowsums[colnum])+1,2)**gamma)))

def resample(prc_cut,Y,Gtmp,train,gamma):
	add = list()
	rem = list()
	colsums = np.sum(Y[train,:],axis=0);
	cutoff = math.ceil(np.percentile(colsums,prc_cut));
	for i in range(len(colsums)):
		if colsums[i] == 0:
			pass
		elif colsums[i] < wt_cutoff(i,cutoff,Gtmp,gamma):
			idx = np.squeeze(np.array(np.where(Y[train,i]>=1)));
			choice = np.random.choice(train[idx],int(wt_cutoff(i,cutoff,Gtmp,gamma)-colsums[i]))
			add = add + choice.tolist();
		elif colsums[i] == wt_cutoff(i,cutoff,Gtmp,gamma):
			pass
		else:
			idx = np.squeeze(np.array(np.where(Y[train,i]>=1)));
			choice = np.random.choice(train[idx],int(colsums[i]-wt_cutoff(i,cutoff,Gtmp,gamma)),replace=False)
			rem = rem + choice.tolist()
	return np.concatenate((list([val for val in train if val not in rem]),add));

def select_feats(Xtmp,num_zero_prc_cut,var_prc_cut):
	#*********************************************************************
	# remove features with many zeros
	num_feat_zeros = np.sum(Xtmp==0,axis=1);
	Xtmp = Xtmp[num_feat_zeros<num_zero_prc_cut*Xtmp.shape[1],:]
	#*********************************************************************
	# remove features with low variance
	feat_vars = np.var(Xtmp,axis=1)
	Xtmp = Xtmp[feat_vars>np.percentile(feat_vars,var_prc_cut),:]
	return(Xtmp)

def get_yn(predict,ys,delta,tau,output_feats):
	D = tf.cast(Dnp, tf.float32);
	G = tf.cast(Gnp, tf.float32);
	ys = tf.cast(ys, tf.float32);
	#print("start")
	Cm = tf.matmul(tf.transpose(tf.matmul(ys,D)),predict+0.1)/tf.reshape(tf.reduce_sum(tf.transpose(tf.matmul(ys,D)),1),(-1,1));
	#print("1")
	mCm = tf.reshape(tf.reduce_mean(tf.cast(tf.matmul(tf.transpose(D),G)>0,tf.float32)*Cm,1),(-1,1));
	#print("2")
	yw = tf.multiply(predict+0.1,tf.matmul(tf.matmul(ys,D),tf.pow(mCm/Cm,tau)));
	#print("3")
	ye = tf.multiply(tf.matmul(ys,G),yw);
	#print("4")
	yt = tf.matmul(ys,tf.matmul(tf.transpose(ys),ye));
	#print("5")
	ya = (delta*yt)+((1-delta)*ye)
	#print("6")
	yn = tf.cast(tf.one_hot(tf.argmax(ya,axis=1),output_feats), dtype=tf.float32)
	#print("7")
	return(yn)

def get_yi(rowsums,G2,ys):
	G2 = tf.cast(G2, tf.float32);
	ys = tf.cast(ys, tf.float32);
	yi = tf.cast(tf.matmul(ys,G2), dtype=tf.float32);
	return(yi)

def run_LAmbDA2(gamma, delta, tau, prc_cut, bs_prc, num_trees, max_nodes):
	global X, Y, Gnp, Dnp, train, test, prt, cv
	D = tf.cast(Dnp, tf.float32);
	G = tf.cast(Gnp, tf.float32);
	#optunity_it = optunity_it+1;
	num_trees = int(num_trees);
	max_nodes = int(max_nodes);
	prc_cut = int(np.ceil(prc_cut));
	print("gamma=%.4f, delta=%.4f, tau=%.4f, prc_cut=%i, bs_prc=%.4f, num_trees=%i, max_nodes=%i" % (gamma, delta, tau, prc_cut, bs_prc, num_trees, max_nodes))
	input_feats = X.shape[1];
	num_labls = G.shape.as_list();
	output_feats = num_labls[1];
	#print(output_feats)
	num_labls = num_labls[0];
	rowsums = np.sum(Gnp,axis=1);
	train2 = resample(prc_cut, Y, Gnp, train, gamma);				# Bug??
	bs = int(np.ceil(bs_prc*train2.size))
	xs = tf.placeholder(tf.float32, [None,input_feats])
	#ys = tf.placeholder(tf.float32, [None,num_labls])
	yin = tf.placeholder(tf.int32, [None])
	print("Vars loaded xs and ys created")
	hparams = tensor_forest.ForestHParams(num_classes=output_feats,
									num_features=input_feats,
									num_trees=num_trees,
									max_nodes=max_nodes).fill()
	print("Tensor forest hparams created")								
	forest_graph = tensor_forest.RandomForestGraphs(hparams)
	print("Tensor forest graph created")
	train_op = forest_graph.training_graph(xs, yin)
	loss_op = forest_graph.training_loss(xs, yin)
	print("Loss and train ops created")
	predict, _, _ = forest_graph.inference_graph(xs)
	print("Tensor forest variables created through predict")
	accuracy_op = tf.reduce_mean(tf.reduce_sum(tf.square(tf.one_hot(yin,output_feats)-predict),reduction_indices=[1]))
	print(tf.reduce_sum(tf.square(tf.one_hot(yin,output_feats)-predict),reduction_indices=[1]))
	#predict = tf.one_hot(pred);
	print("Lambda specific variables created")
	# Creating training and testing steps
	G2 = np.copy(Gnp);
	G2[rowsums>1,:] = 0;
	YI = np.matmul(Y,G2);
	YIrs = np.sum(YI,axis=1);
	trainI = train2[np.in1d(train2,np.where(YIrs==1))];
	print("data type trainI,",trainI.dtype)
	testI = test[np.in1d(test,np.where(YIrs==1))];
	print("trainI testI created")
	#init_vars=tf.global_variables_initializer()
	init_vars = tf.group(tf.global_variables_initializer(),
	resources.initialize_resources(resources.shared_resources()))
	sess = tf.Session()
	sess.run(init_vars)
	print("Session started")
	#beep = sess.run(predict,feed_dict={xs:X[1:100,:]});
	#beep = sess.run(predict,feed_dict={xs:X[train2[0:bs],:]});
	tensor_trainI = {xs: X[trainI, :], yin: sess.run(tf.argmax(get_yi(rowsums,G2,Y[trainI, :]),axis=1))}
	print("tensor_trainI made")
	tensor_testI = {xs: X[testI, :], yin: sess.run(tf.argmax(get_yi(rowsums,G2,Y[testI, :]),axis=1))}
	print("tensor_testI made")
	tensor_train = {xs: X[train2[0:bs], :], yin: sess.run(tf.argmax(get_yn(sess.run(predict,feed_dict={xs:X[train2[0:bs],:]}),Y[train2[0:bs], :],delta,tau,output_feats),axis=1))}
	print("tensor_train made")
	tensor_test = {xs: X[test, :], yin: sess.run(tf.argmax(get_yn(sess.run(predict,feed_dict={xs:X[test,:]}),Y[test, :],delta,tau,output_feats),axis=1))}
	print("tensor_test made")
	#**********************************
	#print("Loss and training steps created with sample tensors")
	# Setting params and initializing
	print("Beginning iterations")
	# Starting training iterations
	print(X.shape)
	for i in range(1,101):
		if i < 50:
			sess.run(train_op, feed_dict=tensor_trainI)
			#print("ran train op")
			if i % 10 == 0:
				print(str(sess.run(accuracy_op, feed_dict=tensor_trainI)) + ' ' + str(sess.run(accuracy_op, feed_dict=tensor_testI)))
		else:
			sess.run(train_op, feed_dict=tensor_train)
			if i % 10 == 0:
				print(str(sess.run(accuracy_op, feed_dict=tensor_train)) + ' ' + str(sess.run(accuracy_op, feed_dict=tensor_test)))
			elif i % 10 == 0:
				np.random_shuffle(train2);
				tensor_train = {xs: X[train2[0:bs], :], yin: sess.run(get_yn(sess.run(predict,feed_dict={xs:X[train2[0:bs],:]}),Y[train2[0:bs], :],delta,tau,output_feats))}
	if prt:
		blah = sess.run(predict, feed_dict=tensor_test);
		sio.savemat('preds_cv' + str(cv) + '.mat', {'preds': blah});
		sio.savemat('truth_cv' + str(cv) + '.mat', {'labels': Y[test, :]});
	acc = sess.run(accuracy_op, feed_dict=tensor_test) 
	print("loss1=%.4f, gamma=%.4f, delta=%.4f, tau=%.4f, prc_cut=%i, bs_prc=%.4f, num_trees=%i, max_nodes=%i" % (acc, gamma, delta, tau, prc_cut, bs_prc, num_trees, max_nodes))
	tf.reset_default_graph();
	return(acc)
