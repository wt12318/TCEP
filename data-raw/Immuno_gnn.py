import torch_geometric
import torch
import pandas as pd 
from pMHCDataset import pMHCDataset
from torch_geometric.loader import DataLoader
from torch_geometric.data import Dataset, Data
import numpy as np
import os
import shutil
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-a", "--aaindex", help="The aaindex file path")
parser.add_argument("-m", "--model", help="The model .pth file path")
parser.add_argument("-t", "--threads", help="The number of threads", type=int)
parser.add_argument("-i", "--input", help="The input file")
parser.add_argument("-o", "--outdir", help="The dir of output file")

args = parser.parse_args()
aaindex_file = os.path.abspath(os.path.expanduser(args.aaindex))##去掉可能的 ~
model_file = os.path.abspath(os.path.expanduser(args.model))
threads_num = args.threads 
out_dir = os.path.abspath(os.path.expanduser(args.outdir))
input_file = os.path.abspath(os.path.expanduser(args.input))
input_file_name = os.path.basename(input_file)

aaindex = pd.read_csv(aaindex_file)
model = torch.load(model_file)

dir_root = out_dir + "/data/"
dir_raw = out_dir + "/data/raw/"
dir_process = out_dir + "/data/processed/"

if (os.path.exists(dir_raw)):
    shutil.rmtree(dir_raw, ignore_errors=True)
os.makedirs(dir_raw)
if (os.path.exists(dir_process)):
    shutil.rmtree(dir_process, ignore_errors=True)
os.makedirs(dir_process)
shutil.copy2(input_file,dir_raw)

torch.set_num_threads(threads_num)
test_dataset = pMHCDataset(root=dir_root,filename=input_file_name,aaindex=aaindex,test=False, val=False)
test_dt = pd.read_csv(input_file)

dt_loader = DataLoader(test_dataset, batch_size=128, shuffle=False)
model.eval()
model = model.double()
all_preds = []
all_labels = []
all_preds_raw = []
for batch in dt_loader:
    pred = model(batch.x.double(),batch.batch)
    all_preds.append(np.rint(torch.sigmoid(pred).cpu().detach().numpy()))
    all_labels.append(batch.y.cpu().detach().numpy())
    all_preds_raw.append(torch.sigmoid(pred).cpu().detach().numpy())
all_preds = np.concatenate(all_preds).ravel()
all_labels = np.concatenate(all_labels).ravel()
all_preds_raw = np.concatenate(all_preds_raw).ravel()
df = pd.DataFrame(zip(all_labels, all_preds_raw, all_preds), columns=['label','pred_raw','pred'])

df_a = pd.concat([test_dt.reset_index(drop=True), df], axis=1)
df_a.to_csv(out_dir + "/pred_res.csv")