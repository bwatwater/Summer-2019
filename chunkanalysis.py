import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

''' Bennett Atwater- Summer 2019
Program for analyzing standard deviation and mean drift over time. Good for analyzing long-term data scans. Take note of sample 
rate and period of time scanned to determine function parameters. 
'''


path=r'C:\Users\Bennett\Desktop\sensitivity3alt.txt'
#n=number of points per chunk. Related to period of time per chunk.
#m=factor to decrease number of points. Good for tuning sample rate. Occurs after splitting up the data into chunks. 
#k=total number of chunks
fig, ax = plt.subplots()
means=[]
stdlist=[]
chunk=[]
def chunks(n,m,k):
    plt.figure(1)
    for i in range(k):
        data=pd.read_csv(path, sep=" ", header=None, names=['t','v'], skiprows=i*n, nrows=n)
        newdata=data.iloc[lambda x: x.index % m==0]
        t=newdata.t
        v=newdata.v
        stdlist.append(np.std(v))
        ax.plot(t,v, 'r.', alpha=0.5, label='chunk '+str(i+1))
        chunk.append(i+1)
        means.append(np.mean(v))
    print(stdlist)
    print(np.mean(stdlist))
    figure(1)
    plt.xlabel('Time (sec)')
    plt.ylabel('Photodiode Signal (Arb)')
    plt.text(-5000, 9, 'mean of 1-minute standard deviations: '+str(np.mean(stdlist)))
    plt.title('Photodiode Signal versus Time (Stability Test 2)')
    figure
    plt.figure(2)
    plt.hist(stdlist)
    plt.xlabel('Standard Deviation (signal units)')
    plt.ylabel('Probability')
    plt.title('Standard Deviation Histogram (Stability Test 2)')
    plt.figure(3)
    plt.plot(chunk,stdlist, 'r.')
    plt.xlabel('Chunk number')
    plt.ylabel('Standard deviation (Signal Units)')
    plt.title('Standard deviation versus chunk')
    figure(4)
    plt.plot(chunk,means, 'r.')
    plt.xlabel('Chunk number')
    plt.ylabel('Mean over chunk period (Signal Units)')
    plt.title('Mean signal value versus chunk')
    plt.show()


chunks(120000,10,70)

