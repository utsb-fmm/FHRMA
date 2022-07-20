#Train FSMHR model

#Associated paper in https://www.preprints.org/manuscript/202207.0131/v1

#FHR Morphological Analysis Toolbox Copyright (C) 2022 Samuel Boudet, Faculté de Médecine et Maïeutique, samuel.boudet@gmail.com
#This file is part of FHR Morphological Analysis Toolbox
#FHR Morphological Analysis Toolbox is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#FHR Morphological Analysis Toolbox is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/.


import tensorflow as tf
import array
import numpy as np
import os
import time
import datetime
import pickle
from IPython.core.debugger import set_trace
from tensorflow.keras import layers as lay
from tensorflow.keras import regularizers
from official.common import distribute_utils


############

Mname='FSScalp'
batch_size=40
grufn=lay.GRU
EpochPerBatch=[8,8]
#8core*2steps*5batchsize

basefolder=""
logfolder="gs://fhrfalsesigus/logs/"

############
tf.profiler.experimental.server.start(9012)
strategy = distribute_utils.get_distribution_strategy(
  distribution_strategy="tpu",
  num_gpus=1,      
  tpu_address="fhrfalsesigpre1")
  
strategy_scope = distribute_utils.get_strategy_scope(strategy)
  
print("REPLICAS: ", strategy.num_replicas_in_sync)


############


with open(basefolder+'dataV8int40-40.pkl','rb') as f:
    X_train, Y_train, X_val, Y_val, ListRec,LengthRec = pickle.load(f)
samp=X_train.shape[1]

Y_train[:,:,0]=Y_train[:,:,0]*EpochPerBatch[0]
Y_val[:,:,0]=Y_val[:,:,0]*EpochPerBatch[1]
print("Data Loaded")
############


@tf.function 
def weighted_binary_crossentropy(yTrue, yPred): 
    Pred=yPred[:,:,0]
    Pred=tf.clip_by_value(Pred,tf.constant(.00001),tf.constant(.999999))
    isFalse=yTrue[:,:,1]
    docare=yTrue[:,:,0]    
    N=(tf.keras.backend.log(Pred)*isFalse+tf.keras.backend.log(tf.constant(1.)-Pred)*(tf.constant(1.)-isFalse))*docare
    return -tf.keras.backend.sum(N)#/tf.keras.backend.sum(docare)

@tf.function 
def weighted_accuracy(yTrue, yPred):
    
    Pred=yPred[:,:,0]
    isFalse=yTrue[:,:,1]
    docare=yTrue[:,:,0] 
    C=tf.cast((Pred>tf.constant(.5)),tf.float32)    
    N=tf.keras.backend.sum( tf.cast((C==isFalse),tf.float32) * docare)
    return tf.keras.backend.sum(N)#/tf.keras.backend.sum(docare)

############


class Split(tf.keras.layers.Layer):
    def __init__(self):
        super(Split, self).__init__()      
    def call(self, inputs):
        return tf.split(inputs,num_or_size_splits=2, axis=0)


class RevT(tf.keras.layers.Layer):
    def __init__(self):
        super(RevT, self).__init__()      
    def call(self, inputs):
        return tf.keras.backend.reverse(inputs,1)  

class Set0forReset(tf.keras.layers.Layer):
    def __init__(self):
        super(Set0forReset, self).__init__()   

    def call(self, inputs):
        return tf.concat([tf.multiply(inputs[:,:,0:-1],1.-inputs[:,:,-1,None]),inputs[:,:,-1,None]],axis=2)        
        
################
def siglossgenerator(X,Y):
    pFHRn=np.array([4,12,24,36,64,120])
    pFHRp=np.array([.25,.25,.25,.25,.25,.25])
    pFHRd=np.array([2000,500,150,60,20,6])
    

    pKeepFHR=.20
    pCut=.3
    pNoStage2=.1
    stdMult=tf.constant(.08,dtype=tf.float32)

   
    rg=tf.reshape(tf.range(0,samp,dtype=tf.int32),[samp,1])
    if tf.random.uniform([])<pCut:
        pos=tf.random.uniform([],maxval=samp-1,dtype=tf.int32)
        newX=tf.tensor_scatter_nd_update(X, [[pos,0],[pos,1],[pos,2],[pos,3]], [0,0,0,1])
        newY=tf.tensor_scatter_nd_update(Y, [[pos,0],[pos,1]], [0,0])
    else:
        newX=X
        newY=Y  

    FHR=newX[:,0:1]
    maskFHR=newX[:,1:2]
    isStage2=newX[:,2:3]
    care=newY[:,0:1]
    falsesig=newY[:,1:2]

    if tf.random.uniform([])<pNoStage2:
        isStage2=tf.zeros(tf.shape(isStage2),dtype=tf.dtypes.float32)  
    
    if tf.random.uniform([])>pKeepFHR:
        NoiseFactor=tf.math.pow(2.,tf.random.normal([],-1,1,dtype=tf.float32))
        for k in range(len(pFHRn)):
            N=tf.cast(tf.floor(NoiseFactor*pFHRn[k]),tf.int32)
            for j in range(N):
                if tf.random.uniform([])<pFHRp[k]: # removing maternal heart rate
                    mid=tf.random.uniform([],maxval=samp,dtype=tf.int32)
                    dur=tf.cast(abs(tf.random.normal([],mean=pFHRd[k],stddev=pFHRd[k]/2)),tf.int32)
                    maskFHR=tf.where( (rg<mid-dur)|(rg>=mid+dur), maskFHR, [0]) 
                    care=tf.where( (rg<mid-dur)|(rg>=mid+dur), care, [0]) 

        N=tf.cast(tf.floor(NoiseFactor*30),tf.int32)
        for j in range(N):
            mid=tf.random.uniform([],maxval=samp,dtype=tf.int32)
            durT=tf.cast(abs(tf.random.normal([],mean=120,stddev=120/2)),tf.int32)
            durPart=tf.cast(abs(tf.random.normal([],mean=60,stddev=60/2)),tf.int32)
            OffsetPart=tf.cast(tf.random.normal([],mean=0,stddev=20),tf.int32)
            endFalse=tf.clip_by_value(durPart+OffsetPart,-durT,durT);
            startFalse=tf.clip_by_value(-durPart+OffsetPart,-durT,durT);
            r=tf.random.uniform([])
            change=tf.where( (rg>=mid+startFalse) & (rg<mid+endFalse) & (falsesig==0.)& (care>0.), [1.], [0.] )
            
            falsesig=falsesig+change
            maskFHR=tf.where( (rg<mid-durT)|(rg>=mid+durT)|(change==1.), maskFHR, [0])
            care=tf.where( (rg<mid-durT)|(rg>=mid+durT)|(change==1.), care*(1+change), [0]) 
            if r<.5:
                FHR=(FHR+2)*(1+change)-2
            else:
                FHR=(FHR+2)*(1-.5*change)-2
  


    maskFHR=tf.where( FHR<=2.25, maskFHR, [0])   
    
    
    S=1.+tf.random.normal([],stddev=stdMult)

    newX=tf.concat([ ((FHR+2)*S-2)*maskFHR , maskFHR,isStage2,newX[:,3:4] ],axis=1)
    
    newX=tf.ensure_shape(newX,X.shape)
    newY=tf.concat([ care,falsesig ],axis=1)
    return newX,newY



traindataset=tf.data.Dataset.from_tensor_slices(    (  tf.constant(X_train,dtype=tf.float32)  ,  tf.constant(Y_train,dtype=tf.float32)  )    )
traindataset=traindataset.shuffle(X_train.shape[0],reshuffle_each_iteration=True).repeat()
traindataset=traindataset.map(siglossgenerator)
traindataset=traindataset.batch(batch_size,drop_remainder=True).prefetch(tf.data.experimental.AUTOTUNE)

valdataset=tf.data.Dataset.from_tensor_slices((X_val,Y_val)).cache().repeat().batch(batch_size,drop_remainder=True)
################
class Sparsity(tf.keras.constraints.Constraint):
    def __init__(self, n,h):
        self. n = n    
        self.h=h
        m=np.zeros([h,h])
        for i in range(int(h/n)):
            m[i*n:i*n+n,i*n:i*n+n]=1.
        self.mask=tf.constant(np.concatenate([m,m,m],axis=1),dtype=tf.float32)
        
    def __call__(self, w):
        w.assign(w*self.mask)
        return w

class SparsityKernel(tf.keras.constraints.Constraint):
    def __init__(self, m,n,slices):
        #m input size ; n state size
        #slices
        ma1=np.zeros([m,n],dtype=np.float32)
        for i in range(slices.shape[0]):
            ma1[slices[i,0]:slices[i,1],slices[i,2]:slices[i,3]]=1.
        self.mask=tf.constant(np.concatenate([ma1,ma1,ma1],axis=1),dtype=tf.float32)
        ma2=np.zeros([m,3*n])
        ma2[m-1,n:2*n]=-1.e30
        ma2[m-1,0:n]=-1.e30
        self.maskReset=tf.constant(ma2,dtype=tf.float32)
    def __call__(self, w):
        w.assign(w*self.mask+self.maskReset)
        return w

class Reseter(tf.keras.constraints.Constraint):
    def __init__(self, m,n):
        ma=np.zeros([m,3*n])
        ma[m-1,n:2*n]=-1.e30
        ma[m-1,0:n]=-1.e30
        self.mask=tf.constant(ma,dtype=tf.float32)
        
    def __call__(self, w):
        w.assign(w+self.mask)
        return w

class BiasConstraintCallback(tf.keras.callbacks.Callback):
    def on_train_batch_end(self, batch, logs=None):
        for GRUlayer in ["GRU1","GRU2","GRU3"]:
            W=self.model.get_layer(GRUlayer).get_weights()
            L=int(len(W[0][-1])/3)
            W[0][-1][2*L:3*L]=-W[2][0][2*L:3*L]#-W[2][1][2*L:3*L] #Pour CudNN a priori W[2][2*L:3*L]-W[2][5*L:6*L] Mais à vérifier ordre des 3 couches
            self.model.get_layer(GRUlayer).set_weights(W)
            
################
                    

def create_model(batch_size,siglen):
    I=lay.Input(batch_shape= (batch_size, siglen, 4), name="SigInput") 

    RevI=RevT()(I)
    IforwBack=lay.concatenate([I,RevI], axis=0)
    
    GRU1Dop=grufn(24,return_sequences=True,recurrent_initializer='glorot_uniform'
        ,stateful=False,recurrent_constraint=Sparsity(2,24),name='GRU1')(IforwBack)
    GRU1Dop=lay.GaussianDropout(.2)(GRU1Dop)
    Lay1Forw,Lay1Backw=Split()(GRU1Dop)
    Lay1=Set0forReset()(     lay.concatenate(  [lay.concatenate([Lay1Forw,RevT()(Lay1Backw),I], axis=2), lay.concatenate([RevT()(Lay1Forw),Lay1Backw,RevI], axis=2)],axis=0  )     )
    

    slicesKernel2=np.array([[0,4,0,6],[12,16,0,6],[24,28,0,6],[36,40,0,6]
        ,[4,8,6,12],[16,20,6,12],[28,32,6,12],[40,44,6,12]
        ,[8,12,12,18],[20,24,12,18],[32,36,12,18],[44,48,12,18]
        ,[48,52,0,18]])
    GRU2Dop=grufn(18,return_sequences=True,recurrent_initializer='glorot_uniform'
        ,stateful=False,recurrent_constraint=Sparsity(3,18),kernel_constraint=SparsityKernel(52,18,slicesKernel2),name='GRU2')(Lay1)
    GRU2Dop=lay.GaussianDropout(.3)(GRU2Dop)
    Lay2Forw,Lay2Backw=Split()(GRU2Dop)
    Lay2=Set0forReset()(     lay.concatenate(  [lay.concatenate([Lay2Forw,RevT()(Lay2Backw),I], axis=2), lay.concatenate([RevT()(Lay2Forw),Lay2Backw,RevI], axis=2)],axis=0  )      )
    
    GRU3Dop=grufn(12,return_sequences=True,recurrent_initializer='glorot_uniform'
        ,stateful=False,recurrent_constraint=Sparsity(4,12),kernel_constraint=Reseter(40,12),name='GRU3')(Lay2)
    GRU3Dop=lay.GaussianDropout(.4)(GRU3Dop)
    Lay3Forw,Lay3Backw=Split()(GRU3Dop)
    CDop=lay.concatenate([Lay3Forw,RevT()(Lay3Backw)], axis=2)
    
    PDop=lay.Dense(1,activation='sigmoid',name="Dense1")(CDop)

    model = tf.keras.Model(inputs=[I], outputs=[PDop])

    return model

################
    
callbacks=[tf.keras.callbacks.TensorBoard(log_dir = logfolder+Mname , histogram_freq = 0)]

callbacks.append(BiasConstraintCallback())

callbacks.append(tf.keras.callbacks.ModelCheckpoint(filepath=Mname+"/best-{epoch:05d}-{val_loss:.4f}.h5",
                                                 save_best_only=True,
                                                 verbose=1))
callbacks.append(tf.keras.callbacks.ModelCheckpoint(filepath=Mname+"/save-{epoch:05d}.h5",
                                                 save_freq=200,
                                                 verbose=1))
def scheduler(epoch):
    if epoch<100:
        return 5.e-3
    elif epoch<300:
        return 2.e-3
    elif epoch<1000:
        return 1.e-3
    else:
        return 5.e-4
callbacks.append(tf.keras.callbacks.LearningRateScheduler(scheduler))


################


initial_epoch=0
with strategy_scope:
    model=create_model(int(batch_size/strategy.num_replicas_in_sync),samp)
    #model.load_weights(basefolder+Mname+"/save-00800.h5")
    model.compile(
        optimizer=tf.optimizers.Adam(learning_rate=1e-3),
        loss=weighted_binary_crossentropy,
        metrics=[weighted_accuracy])
        
        
model.fit(traindataset,
    epochs=20000,
    steps_per_epoch=np.floor(X_train.shape[0]/batch_size),
    validation_freq=1,
    initial_epoch=initial_epoch,
    validation_steps=1,
    validation_data=valdataset,
    callbacks=callbacks    
)
