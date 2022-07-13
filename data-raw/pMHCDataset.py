# -*- coding: utf-8 -*-
import pandas as pd
import torch
import torch_geometric
from torch_geometric.data import Dataset, Data
import numpy as np 
import os
from tqdm import tqdm
import copy
class pMHCDataset(Dataset):
    def __init__(self, root, filename, aaindex, test=False, val=False, transform=None, pre_transform=None):
        """
        root = Where the dataset should be stored. This folder is split
        into raw_dir (downloaded dataset) and processed_dir (processed data). 
        """
        self.test = test
        self.val = val
        self.filename = filename
        #self.pep_tcr_mat = pep_tcr_mat_dir
        #self.hla_tcr_mat = hla_tcr_mat_dir
        self.aaindex = aaindex
        #self.aa_code = aa_code
        super(pMHCDataset, self).__init__(root, transform, pre_transform)
        
    @property
    def raw_file_names(self):
        """ If this file exists in raw_dir, the download is not triggered.
            (The download func. is not implemented here)  
        """
        return self.filename

    @property
    def processed_file_names(self):
        """ If these files are found in processed_dir, processing is skipped"""
        self.data = pd.read_csv(self.raw_paths[0]).reset_index()
        
        if self.test:
            return [f'data_test_{i}.pt' for i in list(self.data.index)]
        elif self.val:
            return [f'data_val_{i}.pt' for i in list(self.data.index)]
        else:
            return [f'data_{i}.pt' for i in list(self.data.index)]

    def download(self):
        pass##不需要下载
    
    ##获取后部分的 edge_index
    #def _get_index(self, raw_edge):
    #    raw_edge_copy = copy.deepcopy(raw_edge)
    #    for j in raw_edge_copy[1]:
    #        if j in raw_edge[0]:
    #            continue
    #        b = list(compress(raw_edge_copy[0], [int(i == j) for i in raw_edge_copy[1]]))##哪些是j的source 节点
    #        raw_edge[0].extend([j]*len(b))
    #        raw_edge[1].extend(b)
    #    return raw_edge
    
    def process(self):
        self.data = pd.read_csv(self.raw_paths[0])
        for index, sample in tqdm(self.data.iterrows(), total=self.data.shape[0]):#tqdm可以显示运行进程
            # Get node features
            node_feats = self._get_node_features(sample["pep"],sample["hla_seq"],self.aaindex)
            # Get adjacency info
            #edge_index = self._get_adjacency_info(sample["pep"],sample["hla_seq"],self.pep_hla_mat)
            ##Get edge weight
            #edge_weight = self._get_edge_weight(sample["pep"],edge_index,self.pep_hla_mat)
            # Get labels info
            label = self._get_labels(sample["type"])

            # Create data object
            data = Data(x=node_feats, 
                        #edge_index=edge_index,
                        y=label
                        #edge_weight=edge_weight
                        ) 
            if self.test:
                torch.save(data, 
                    os.path.join(self.processed_dir, 
                                 f'data_test_{index}.pt'))
            elif self.val:
                torch.save(data, 
                    os.path.join(self.processed_dir, 
                                 f'data_val_{index}.pt'))                
            else:
                torch.save(data, 
                    os.path.join(self.processed_dir, 
                                 f'data_{index}.pt'))

    def _get_node_features(self,pep,HLA,aaindex):
        """ 
        This will return a matrix / 2d array of the shape
        [Number of Nodes, Node Feature size]
        """
        all_seq = pep + HLA
        all_node_feats = []
        for index, aa in enumerate(all_seq):
            node_feats = []
            ##aaindex
            node_feats.extend(aaindex[aa].to_list())
            ##氨基酸的 one-hot 编码
            #onehot = [0] * len(aa_code)
            #onehot[aa_code.index(aa)] = 1
            #node_feats.extend(onehot)
            ##氨基酸属于哪个序列
            anchar = [0,len(pep)]
            seq_onehot = [0,0]
            seq_onehot[sum([index >= i for i in anchar])-1] = 1
            node_feats.extend(seq_onehot)
            all_node_feats.append(node_feats)
        all_node_feats = np.asarray(all_node_feats)
        #min_max_scaler = preprocessing.MinMaxScaler()
        #all_node_feats_trans = min_max_scaler.fit_transform(all_node_feats)##min-max 转化
        return torch.tensor(all_node_feats, dtype=torch.float)

    def _get_adjacency_info(self,pep,HLA,pep_hla):
        """
        We could also use rdmolops.GetAdjacencyMatrix(mol)
        but we want to be sure that the order of the indices
        matches the order of the edge features
        """
        ##生成边
        edge_index = [[],[]]
        pep_hla_mat = pd.read_csv(pep_hla+"mat_"+HLA+".csv")
        for index, node in enumerate(pep):
            b = [i > 0.7 for i in pep_hla_mat.loc[:,"p"+str(index)].tolist()]
            which_node = [i+len(pep) for i, x in enumerate(b) if x]
            edge_index[0].extend([index]*len(which_node))
            edge_index[1].extend(which_node)
        edge_index = self._get_index(edge_index)

        for index, node in enumerate(pep):
            if index == 0 or index == len(pep)-1:##节点是第一个或最后一个
                edge_index[0].extend([index]*1)
                if index == 0:##第一个和后面一个有连接
                    edge_index[1].extend([index+1])
                else:##最后一个和前面一个有连接
                    edge_index[1].extend([index-1])
            else:##中间的和前面后面都有连接
                edge_index[0].extend([index]*2)
                edge_index[1].extend([index-1,index+1])
        for index, node in enumerate(HLA):
            if index == 0 or index == len(HLA)-1:##节点是第一个或最后一个
                edge_index[0].extend([index+len(pep)]*1)
                if index == 0:##第一个和后面一个有连接
                    edge_index[1].extend([index+len(pep)+1])
                else:##最后一个和前面一个有连接
                    edge_index[1].extend([index+len(pep)-1])
            else:##中间的和前面后面都有连接
                edge_index[0].extend([index+len(pep)]*2)
                edge_index[1].extend([index+len(pep)-1,index+len(pep)+1])
        edge_index = torch.tensor(edge_index)
        return edge_index
        
    def _get_labels(self, label):
        label = np.asarray([label])
        return torch.tensor(label, dtype=torch.int64)

    def len(self):
        return self.data.shape[0]

    def get(self, idx):
        """ - Equivalent to __getitem__ in pytorch
            - Is not needed for PyG's InMemoryDataset
        """
        if self.test:
            data = torch.load(os.path.join(self.processed_dir, 
                                 f'data_test_{idx}.pt'))
        elif self.val:
            data = torch.load(os.path.join(self.processed_dir, 
                                 f'data_val_{idx}.pt'))            
        else:
            data = torch.load(os.path.join(self.processed_dir, 
                                 f'data_{idx}.pt'))   
        return data
