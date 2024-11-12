import sys
import pandas as pd
from multiprocessing import Pool, cpu_count

sys.path.append('../be_predict_efficiency')
sys.path.append('../be_predict_bystander')

from be_predict_efficiency import predict as be_efficiency_model
from be_predict_bystander import predict as be_bystander_model


#From pandas, import annotated_csv
guideDf = pd.read_csv("./annotated_guides.csv")


#Load model
be_efficiency_model.init_model(base_editor = "ABE", celltype = "mES")
be_bystander_model.init_model(base_editor="CBE", celltype="mES")


'''
Sample usage:
seq = 'TATCAGCGGGAATTCAAGCGCACCAGCCAGAGGTGTACCGTGGACGTGAG'
pred_d = be_efficiency_model.predict(seq)

'''

#Potentially multiprocessing to speed things up? Unsure if there is a batch one, check with calvin
def predict_wrapper(seq):
    return be_efficiency_model.predict(seq)

def byastander_wrapper(seq):
    return be_efficiency_model.predict(seq)


guides = guideDf['sgRNA_seq'].tolist()
#Creation of thread pool, mapping predict to all guides list
pool = Pool(cpu_count())
results = pool.map(predict_wrapper, guides)
pool.close()
pool.join()
