from tqdm import tqdm
from sed import SEDFit
import numpy as np
import tensorflow as tf
from sklearn.metrics import confusion_matrix

x=SEDFit(grid_type='btsettl')
l=len(x.grid_table)
w=len(x.sed)
la=np.log10(x.sed['la'])

def makedata():
    arr=np.zeros((l,w))
    for i in tqdm(range(l)):
        x.addguesses(r=[1],teff=x.grid_table['teff'][i],logg=x.grid_table['logg'][i],feh=x.grid_table['feh'][i],av=np.random.beta(0.3,5)*20)
        arr[i]=x.mags
    return arr


def reshape(data,label):
    xx=np.tile(data,(w,1))
    y=np.array([np.tile(np.arange(w),(l,1)).T.flatten()]).T
    yy=np.concatenate([xx,y],axis=1)
    
    zz=label.T.flatten()
    
    return yy,zz

def worsten(data):
    noise=np.random.lognormal(mean=-4, sigma=1,size=(l,))
    noise=np.transpose(np.tile(noise,(w,1)))
    gauss=np.random.normal(size=(l, w))
    data=data+gauss*noise
    for i in range(l):
        random=np.round(np.random.beta(np.random.uniform(0.01,4.3),1,size=w))
        data[i]=data[i]*random
        a=np.where(data[i]!=0)[0]
        data[i,a]=data[i,a]-np.mean(data[i])
        
    label=np.zeros_like(data)+2
    a=np.where(data[i]==0)
    label[a]=0

    for i in range(l):
        a=np.where(data[i]!=0)[0]
        if np.random.uniform()>0.5:
            bad=int(np.random.uniform(low=0,high=np.min([5,len(a)])))
            badid=np.random.choice(a,size=bad,replace=False)
            offset=np.random.lognormal(0.01,size=bad)*np.sign(np.random.uniform(low=-0.1,high=0.1,size=bad))
            data[i,badid]=data[i,badid]+offset
            label[i,badid]=1
            
    data,label=reshape(data,label)
    return data,label
    
data=makedata()

xx,zz=[],[]
for i in range(200):
    da,lb=worsten(data.copy())
    xx.append(da)
    zz.append(lb)
xx=np.concatenate(xx)
zz=np.concatenate(zz)
xx,zz=reshape(xx,zz)

xx1,zz1=worsten(data.copy())
xx1,zz1=reshape(xx1,zz1)

a=np.where(zz==0)[0]
b=np.where(zz==1)[0]
c=np.where(zz==2)[0]
d=np.where(zz==3)[0]
a=np.random.choice(a,len(b),replace=False)
d=np.random.choice(d,len(b),replace=False)
e=np.concatenate([a,b,c,d])
e=np.random.choice(e,len(e),replace=False)
xx,zz=xx[e],zz[e]



tf.keras.backend.clear_session()
input1= tf.keras.Input(shape=(w,2))

layers1 = tf.keras.layers.Conv1D(16, 2, activation='tanh',padding="same")(input1)
layers1 = tf.keras.layers.MaxPool1D(2)(layers1)
layers1 = tf.keras.layers.Conv1D(32, 2, activation='tanh',padding="same")(layers1)
layers1 = tf.keras.layers.MaxPool1D(2)(layers1)
layers1 = tf.keras.layers.Conv1D(32, 2, activation='tanh',padding="same")(layers1)
layers1 = tf.keras.layers.MaxPool1D(2)(layers1)
layers1 = tf.keras.layers.Conv1D(64, 2, activation='tanh',padding="same")(layers1)
layers1 = tf.keras.layers.MaxPool1D(2)(layers1)
layers1 = tf.keras.layers.Flatten()(layers1)
layers1 = tf.keras.layers.Dense(100,activation='tanh')(layers1)
layers1 = tf.keras.layers.Dense(100,activation='tanh')(layers1)
layers1 = tf.keras.layers.Dense(100,activation='tanh')(layers1)
output = tf.keras.layers.Dense(4,activation='softmax')(layers1)

model = tf.keras.Model(inputs=input1, outputs=output)
model.compile(optimizer=tf.keras.optimizers.Adam(learning_rate=1e-4),
              loss=tf.keras.losses.SparseCategoricalCrossentropy(from_logits=True),
              metrics=['accuracy'])
model.summary()

model.fit(xx,zz, validation_data=(xx1,zz1), epochs=15000,batch_size=50000,
          callbacks=[tf.keras.callbacks.EarlyStopping(monitor='val_accuracy',patience=10)])
          
model.save('quality.keras')


predictions = tf.argmax(model.predict(xx1,batch_size=10000),axis=1).numpy()
accuracy = confusion_matrix(zz1, predictions)
print(accuracy)